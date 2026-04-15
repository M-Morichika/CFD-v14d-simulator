// config/experiments.js
export const params = {
    // 空間・時間解像度
    gridSize: 60,      // 100から60へ（計算負荷と安定性のため）
    dt: 0.002,         // 0.05から0.002へ（拡散項の爆発を防ぐ）
    
    // 密度場 ρ のパラメータ (Section 3.1)
    lambda0: 1.0,
    kappa: 0.5,
    D0: 0.05,
    pi0: 1.0,
    
    // 内部代謝 π_int(t) の動的更新 (Section 5.1)
    pi_base: 2.0,
    pi_tgt: 0.2,
    k_rec: 0.05,
    k_fb: 0.3,
    
    // 境界条件 (Section 3.2)
    pi_ext: 1.0,
    pi_act: 0.1,
    
    // 防衛的粘性 γ(G) (Section 3.3)
    G_th: 2.0,
    a_steep: 5.0,
    gamma_min: 1.0,
    gamma_max: 50.0,
    B: 1.0,
    
    // ポテンシャル F の更新 (Section 3.5)
    tau_base: 100.0,
    r0: 0.05,
    rho0: 1.0,
    k_F: 2.0
};