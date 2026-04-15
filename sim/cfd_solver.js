// sim/cfd_solver.js
// CFD v14d: 統合ソルバー (Advection-Diffusion-Reaction + 動的ポテンシャル + 内部代謝の非自律系更新)

import { calc_reaction } from '../core/reaction.js';
import { calc_diffusion } from '../core/diffusion.js';
import { calc_advection } from '../core/advection.js';
import { calc_noise } from '../core/noise.js';
import { calc_G, calc_gamma } from '../core/resistance.js';
import { calc_erosion } from '../core/erosion.js';

export class CFDSolver {
    constructor(config) {
        this.config = config;
        this.N = config.gridSize || config.N || 100;
        
        // 状態変数の初期化
        this.rho = new Array(this.N).fill(0).map(() => new Array(this.N).fill(0));
        this.F = new Array(this.N).fill(0).map(() => new Array(this.N).fill(0));
        this.G = 0;
        this.gamma = config.gamma_min || 1.0;

        // v14d: 内部代謝 pi_int を動的変数として管理
        this.pi_int = config.pi_base || 2.0;

        const dx = config.dx || (2.0 / this.N);
        this.dx = dx;
        const half = this.N / 2;
        
        const k_F = config.k_F || 2.0;
        for (let i = 0; i < this.N; i++) {
            for (let j = 0; j < this.N; j++) {
                const x = (j - half + 0.5) * dx;
                const y = (i - half + 0.5) * dx;
                // 進化で固定された長期ポテンシャル F0
                this.F[i][j] = 0.5 * k_F * (x * x + y * y);
            }
        }

        // 浸食用の W-kernel (メキシカンハット型) の事前計算
        const W_sigma = config.W_sigma || 0.1;
        this.k_rad = Math.ceil(W_sigma * 2 / dx);
        this.W_kernel = [];
        const two_sig_sq = 2.0 * W_sigma * W_sigma;
        for (let i = -this.k_rad; i <= this.k_rad; i++) {
            let row = [];
            for (let j = -this.k_rad; j <= this.k_rad; j++) {
                let r2 = (i * dx) * (i * dx) + (j * dx) * (j * dx);
                let w = (1.0 - r2 / two_sig_sq) * Math.exp(-r2 / two_sig_sq);
                row.push(w);
            }
            this.W_kernel.push(row);
        }
    }

    // k_drug は外部摂動 (app.js のボタン等から渡される)
    step(dt, k_drug = 0.0) {
        const c = this.config;
        const N = this.N;
        const dx = this.dx;

        // 1. 大域的観測量 M と G の計算
        let M = 0;
        for (let i = 0; i < N; i++) {
            for (let j = 0; j < N; j++) {
                M += this.rho[i][j] * dx * dx;
            }
        }
        this.G = calc_G(this.rho, this.F, N, dx, c.k_G || 1.0);

        // 2. 内部代謝 pi_int(t) の動的更新 (式7: 雪崩効果と這い上がり遅延)
        const f_M = (c.k_fb || 0.3) * M;
        const dpi = -k_drug * (this.pi_int - (c.pi_tgt || 0.2))
                  + (c.k_rec || 0.05) * ((c.pi_base || 2.0) - this.pi_int)
                  + f_M;
        this.pi_int += dpi * dt;

        // 3. 防衛的粘性の更新 (式4)
        this.gamma = calc_gamma(this.G, c.gamma_min || 1.0, c.gamma_max || 50.0, c.a_steep || 5.0, c.G_th || 2.0);
        const v_coeff = (c.B || 1.0) / this.gamma;

        let next_rho = new Array(N).fill(0).map(() => new Array(N).fill(0));
        let next_F = new Array(N).fill(0).map(() => new Array(N).fill(0));
        const half = N / 2;

        // 4. 空間更新ループ (PDEの計算)
        for (let i = 0; i < N; i++) {
            for (let j = 0; j < N; j++) {
                const x = (j - half + 0.5) * dx;
                const y = (i - half + 0.5) * dx;
                const r = Math.sqrt(x * x + y * y);

                // マルコフブランケット(r=1)の外部は計算しない
                if (r > 1.0) {
                    next_rho[i][j] = 0;
                    next_F[i][j] = this.F[i][j];
                    continue;
                }

                const rho_c = this.rho[i][j];
                const rho_u = (i > 0) ? this.rho[i - 1][j] : rho_c;
                const rho_d = (i < N - 1) ? this.rho[i + 1][j] : rho_c;
                const rho_l = (j > 0) ? this.rho[i][j - 1] : rho_c;
                const rho_r = (j < N - 1) ? this.rho[i][j + 1] : rho_c;

                const F_c = this.F[i][j];
                const F_u = (i > 0) ? this.F[i - 1][j] : F_c;
                const F_d = (i < N - 1) ? this.F[i + 1][j] : F_c;
                const F_l = (j > 0) ? this.F[i][j - 1] : F_c;
                const F_r = (j < N - 1) ? this.F[i][j + 1] : F_c;

                // 各PDE項の計算
                const adv = calc_advection(rho_c, rho_u, rho_d, rho_l, rho_r, F_c, F_l, F_r, F_u, F_d, dx, v_coeff);
                const diff = calc_diffusion(rho_c, rho_u, rho_d, rho_l, rho_r, dx, c.D0 || 0.05, c.pi_ext || 1.0, this.pi_int, c.pi0 || 1.0);
                const react = calc_reaction(rho_c, this.pi_int, c.lambda0 || 1.0, c.kappa || 0.5);
                const noise = calc_noise(c.noise_sigma || 0.001, dt);

                // Robin型境界条件 (式2)
                let boundary_flux = 0;
                if (r >= 0.95 && r <= 1.0) {
                    const I_theta = 1.0; // 本モデルでは一様刺激を仮定
                    boundary_flux = (c.pi_ext || 1.0) * I_theta - (c.pi_act || 0.1) * rho_c;
                }

                next_rho[i][j] = Math.max(0, rho_c + (diff + adv + react + boundary_flux) * dt + noise);

                // ポテンシャル F の更新 (式6)
                const F0 = 0.5 * (c.k_F || 2.0) * r * r;
                const r0 = c.r0 || 0.05;
                const rho0 = c.rho0 || 1.0;
                
                // 正則化された有効時定数
                const tau_F = (c.tau_base || 100.0) / (((r + r0) / 1.0) * (1.0 + rho_c / rho0));
                
                const weathering = -(F_c - F0);
                const erosion = calc_erosion(this.rho, j, i, N, this.k_rad, this.W_kernel, c.erosion_coeff || 0.05);

                next_F[i][j] = F_c + ((weathering + erosion) / tau_F) * dt;
            }
        }

        this.rho = next_rho;
        this.F = next_F;
    }
}