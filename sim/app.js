// sim/app.js
// CFD v14d simulator - main application
import { CFDSolver } from './cfd_solver.js';
import { Renderer } from '../viz/render.js';
// 臨床用語を排除した新しい experiments.js (または params.js) から読み込む想定
import { params as defaultConfig } from '../config/experiments.js'; 

let solver, renderer;
let isRunning = true;
let currentConfig;
let isPerturbed = false; // 外部摂動 k_drug のフラグ

let history_M = new Array(5000).fill(0);
let history_pi = new Array(5000).fill(0);
let history_L = new Array(5000).fill(0);
let betaTrail = [];

function loadExperiment() {
    currentConfig = Object.assign({}, defaultConfig);
    solver = new CFDSolver(currentConfig);
    renderer = new Renderer('heatCanvas', 'chartCanvas', 'loopCanvas', currentConfig.N || currentConfig.gridSize);

    history_M.fill(0);
    history_pi.fill(currentConfig.pi_base); // 初期状態は pi_base
    history_L.fill(0);
    betaTrail = [];

    // 初期密度のシード（相転移を観察するための初期ピーク）
    const center = Math.floor(solver.N / 2);
    solver.rho[center][center] = 2.0;
    solver.rho[center - 1][center] = 1.0;
    solver.rho[center + 1][center] = 1.0;
    solver.rho[center][center - 1] = 1.0;
    solver.rho[center][center + 1] = 1.0;

    // パラメータ表示の更新（UI要素がある場合）
    const paramsHTML = Object.entries(currentConfig).map(([k, v]) => `<b>${k}</b>: ${v}`).join(' | ');
    if(document.getElementById('paramDisplay')) {
        document.getElementById('paramDisplay').innerHTML = paramsHTML;
    }
}

function init() {
    // UIイベントのバインド (HTML側に 'btnPerturb' というIDのボタンがある想定)
    const btnPerturb = document.getElementById('btnPerturb');
    if(btnPerturb) {
        // マウスを押している間だけ摂動 (k_drug) を印加する
        btnPerturb.addEventListener('mousedown', () => isPerturbed = true);
        btnPerturb.addEventListener('mouseup', () => isPerturbed = false);
        btnPerturb.addEventListener('mouseleave', () => isPerturbed = false);
    }

    const btnRun = document.getElementById('btnRun');
    if(btnRun) btnRun.addEventListener('click', () => isRunning = !isRunning);
    
    const btnReset = document.getElementById('btnReset');
    if(btnReset) btnReset.addEventListener('click', () => loadExperiment());

    loadExperiment();
    requestAnimationFrame(loop);
}

function loop(timestamp) {
    if (isRunning) {
        // v14d: 外部摂動の決定（UIからの入力）
        // 摂動時は k_drug = 0.5 (最大レート)、それ以外は 0.0
        let current_k_drug = isPerturbed ? 0.5 : 0.0;

        // ソルバーの1ステップ更新（k_drug を渡し、内部で pi_int や M が更新される）
        solver.step(currentConfig.dt, current_k_drug);

        // ソルバーから更新された状態を取得
        let current_pi = solver.pi_int;
        let total_M = 0;
        let max_rho = 0;
        let max_i = solver.N / 2, max_j = solver.N / 2;

        // 密度場 ρ から最大値と総情報量 M を計算
        for (let i = 0; i < solver.N; i++) {
            for (let j = 0; j < solver.N; j++) {
                let val = solver.rho[i][j];
                if (val > max_rho) {
                    max_rho = val;
                    max_i = i; max_j = j;
                }
                total_M += isNaN(val) ? 0 : val;
            }
        }

        // フォーカス局在度 L (旧版の名残として維持)
        let L_val = (total_M > 0.0001) ? (max_rho / total_M) : 0;

        // 履歴の更新
        history_M.push(total_M); history_M.shift();
        history_pi.push(current_pi); history_pi.shift();
        history_L.push(L_val); history_L.shift();

        // 軌跡の保存
        if (max_rho > 0.05) {
            betaTrail.push({ i: max_i, j: max_j });
            if (betaTrail.length > 80) betaTrail.shift();
        } else if (betaTrail.length > 0) {
            betaTrail.shift();
        }

        // レンダリング
        renderer.drawHeatmap(solver.rho, solver.F, max_rho, betaTrail);
        renderer.drawChart(history_M, history_pi, history_L);
        renderer.drawLoop(history_pi, history_M);

        // UI表示の更新
        let piValElement = document.getElementById('pi_val');
        if(piValElement) {
            piValElement.innerText = current_pi.toFixed(3);
            if (current_pi < currentConfig.lambda0) {
                // 点火閾値を下回った場合（雪崩効果発生中）の警告表示
                piValElement.style.color = "#ff5555";
                piValElement.style.fontWeight = "bold";
            } else {
                piValElement.style.color = "#ffffff";
                piValElement.style.fontWeight = "normal";
            }
        }

        let coordDisplay = document.getElementById('coordDisplay');
        if(coordDisplay) {
            if (max_rho > 0.05) {
                let phys_x = (max_j - solver.N / 2 + 0.5) * currentConfig.dx;
                let phys_y = (max_i - solver.N / 2 + 0.5) * currentConfig.dx;
                coordDisplay.innerText = `Focus (x: ${phys_x.toFixed(2)}, y: ${phys_y.toFixed(2)})`;
            } else {
                coordDisplay.innerText = `Focus (none)`;
            }
        }

        let rhoDisplay = document.getElementById('rhoDisplay');
        if(rhoDisplay) {
            rhoDisplay.innerText = 
                `Max ρ: ${max_rho.toFixed(3)} | M: ${total_M.toFixed(3)} | G: ${solver.G.toFixed(3)} | γ: ${solver.gamma.toFixed(3)}`;
        }
    }
    requestAnimationFrame(loop);
}

init();