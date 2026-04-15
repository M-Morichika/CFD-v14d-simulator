// core/advection.js
// CFD v14d: 移流項の計算
// 支配方程式 (1) 第1項: -∇·(ρv)
// 支配方程式 (3) 流速: v = -(B/γ(G))∇F
// 数値的安定性のために風上差分法（Upwind scheme）を採用

export function calc_advection(rho_c, rho_u, rho_d, rho_l, rho_r, F_c, F_l, F_r, F_u, F_d, dx, v_coeff) {
    // v_coeff は B / γ(G) を事前計算して渡す
    const v_max = 5.0; // クーラン条件（CFL条件）違反を防ぐための流速シーリング

    // 1. 各セル境界におけるポテンシャル勾配から流速ベクトルを計算
    let vx_L = Math.max(-v_max, Math.min(v_max, -v_coeff * (F_c - F_l) / dx));
    let vx_R = Math.max(-v_max, Math.min(v_max, -v_coeff * (F_r - F_c) / dx));
    let vy_U = Math.max(-v_max, Math.min(v_max, -v_coeff * (F_c - F_u) / dx));
    let vy_D = Math.max(-v_max, Math.min(v_max, -v_coeff * (F_d - F_c) / dx));

    // 2. 風上差分法によるフラックスの計算（流れの上流側の密度を採用）
    let flux_x_L = (vx_L > 0) ? vx_L * rho_l : vx_L * rho_c;
    let flux_x_R = (vx_R > 0) ? vx_R * rho_c : vx_R * rho_r;
    let flux_y_U = (vy_U > 0) ? vy_U * rho_u : vy_U * rho_c;
    let flux_y_D = (vy_D > 0) ? vy_D * rho_c : vy_D * rho_d;

    // 3. 移流による発散成分を返す
    return (flux_x_L - flux_x_R + flux_y_U - flux_y_D) / dx;
}