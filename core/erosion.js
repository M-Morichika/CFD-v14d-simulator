// core/erosion.js
// CFD v14d: ポテンシャル地形の動的浸食（空間畳み込み）
// 支配方程式 (6) 右辺第2項: 浸食項
// erosion = - c_ero * ∫∫ W(β, β') ρ(β', t) dβ'
// メキシカンハット型核関数 W による密度場 ρ の空間畳み込みを計算する

export function calc_erosion(rho_grid, cx, cy, N, k_rad, W_kernel, c_ero) {
    let sum = 0;
    
    // 畳み込み演算の境界処理（グリッド外縁のはみ出しを防ぐ）
    const min_i = Math.max(0, cy - k_rad);
    const max_i = Math.min(N - 1, cy + k_rad);
    const min_j = Math.max(0, cx - k_rad);
    const max_j = Math.min(N - 1, cx + k_rad);

    // 核関数（W_kernel）を通じた局所的な密度の影響を積算
    for (let i = min_i; i <= max_i; i++) {
        const ky = i - cy + k_rad; // カーネル配列内のYインデックス
        for (let j = min_j; j <= max_j; j++) {
            const kx = j - cx + k_rad; // カーネル配列内のXインデックス
            sum += rho_grid[i][j] * W_kernel[ky][kx];
        }
    }
    
    // 浸食係数（c_ero）を乗じて、地形が削られる方向（負）として返す
    return -c_ero * sum;
}