// core/reaction.js
// CFD v14d: 反応項の計算
// 支配方程式 (1) 第3・4項: (pi_int - lambda_0)*rho - kappa*rho^2
// トランスクリティカル分岐と自己抑制による容量限界

export function calc_reaction(rho_c, pi_int, lambda_0, kappa) {
    // 自己触媒的増殖 (pi_int - lambda_0)*rho と、容量限界による自己抑制 -kappa*rho^2
    return (pi_int - lambda_0) * rho_c - kappa * rho_c * rho_c;
}