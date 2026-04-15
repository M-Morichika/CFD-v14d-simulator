// core/diffusion.js
// CFD v14d: 拡散項の計算
// 支配方程式 (1) 第2項: D_mask * ∇²ρ
// 有効拡散係数 D_mask = D0 / (1 + (pi_ext + pi_int) / pi_0)
// 直交グリッド上の標準的な5点ラプラシアンによる空間2階微分の近似

export function calc_diffusion(rho_c, rho_u, rho_d, rho_l, rho_r, dx, D0, pi_ext, pi_int, pi_0 = 1.0) {
    // 1. 精度パラメータ（外部・内部）に基づく有効拡散係数の計算
    const D_mask = D0 / (1.0 + (pi_ext + pi_int) / pi_0);
    
    // 2. 5点ステンシルによるラプラシアン（空間2階微分）の計算
    const laplacian = (rho_u + rho_d + rho_l + rho_r - 4.0 * rho_c) / (dx * dx);
    
    // 3. 拡散による情報密度の変化量を返す
    return D_mask * laplacian;
}