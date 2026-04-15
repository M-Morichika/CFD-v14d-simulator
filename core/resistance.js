// core/resistance.js
// CFD v14d: 変容抵抗 G と 防衛的粘性 γ(G) の計算
// 支配方程式 (4) 防衛的粘性: γ(G) = γ_min + (γ_max - γ_min) / (1 + exp(-a*(G - G_th)))
// 支配方程式 (5) 変容抵抗: G[ρ,F] = k_G * ∫∫ ρ * F dx dy

// 1. 変容抵抗 G[ρ,F] の計算（直交グリッド上での面積分近似）
export function calc_G(rho, F, N, dx, k_G = 1.0) {
    let sum = 0;
    const half = N / 2;
    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
            const x = (j - half + 0.5) * dx;
            const y = (i - half + 0.5) * dx;
            const r = Math.sqrt(x * x + y * y);
            
            // 有界円盤領域 (r <= 1.0) の内部のみを積分（中心特異点付近は除外）
            if (r > 1.0 || r < 0.01) continue;
            
            sum += rho[i][j] * F[i][j] * dx * dx;
        }
    }
    return k_G * sum;
}

// 2. 防衛的粘性 γ(G) の計算 (シグモイド関数による相転移的スイッチ)
export function calc_gamma(G, gamma_min, gamma_max, a_steep, G_th) {
    const sigmoid = 1.0 / (1.0 + Math.exp(-a_steep * (G - G_th)));
    return gamma_min + (gamma_max - gamma_min) * sigmoid;
}