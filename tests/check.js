#!/usr/bin/env node
// Pruebas de química y termodinámica para molecular-simulation (sin navegador ni dependencias).
"use strict";

const M = require("../js/combustion.js");

let fail = 0;
const ok = (cond, msg) => { console.log((cond ? "✓ " : "✗ ") + msg); if (!cond) fail++; };
const near = (a, b, tol) => Math.abs(a - b) <= tol;

console.log("\n── 1. Termoquímica del metano (Ley de Hess y energías de enlace) ──");
{
  const dHgas = M.hessReactionEnthalpy(false);
  const dHliq = M.hessReactionEnthalpy(true);
  ok(near(dHgas, -802.5, 0.2), `ΔH° combustión CH₄(g) → H₂O(g) = ${dHgas.toFixed(1)} kJ/mol (NIST ≈ −802.3)`);
  ok(near(dHliq, -890.5, 0.2), `ΔH° combustión CH₄(g) → H₂O(l) = ${dHliq.toFixed(1)} kJ/mol (NIST ≈ −890.4)`);

  const bb = M.bondBalance();
  ok(bb.broken === 2648 && bb.formed === 3450 && bb.net === -802,
    `balance de enlaces: +${bb.broken} − ${bb.formed} = ${bb.net} kJ/mol (coherente con Hess)`);
}

console.log("\n── 2. Conservación exacta de átomos (C, H, O) y barreras Eₐ ≥ ΔH° en las 30 reacciones ──");
{
  const all = [...M.REACTIONS, M.UNIMOL_CHO];
  let allBal = true, allFinite = true, allBarrierOk = true;
  for (const rx of all) {
    const sum = (arr, k) => arr.reduce((acc, s) => acc + M.SPECIES[s][k], 0);
    if (sum(rx.p, "C") !== sum(rx.r, "C") || sum(rx.p, "H") !== sum(rx.r, "H") || sum(rx.p, "O") !== sum(rx.r, "O")) {
      allBal = false;
      console.error("Desbalance en:", rx.id);
    }
    if (!Number.isFinite(rx.dH) || rx.Ea < 0 || !Number.isFinite(rx.P0) || rx.P0 <= 0) allFinite = false;
    if (rx.dH > 0 && rx.Ea < rx.dH) {
      allBarrierOk = false;
      console.error("Barrera menor que ΔH en:", rx.id, rx.Ea, "<", rx.dH);
    }
  }
  ok(all.length === 30 && allBal, `las ${all.length} reacciones elementales conservan exactamente C, H y O`);
  ok(allFinite, "todas las reacciones tienen ΔH° finito, Eₐ ≥ 0 y factor P₀ > 0");
  ok(allBarrierOk, "todas las reacciones endotérmicas cumplen Eₐ ≥ ΔH° (sin violaciones termodinámicas)");

  // Verificar valores específicos de reacciones emblemáticas
  const init = M.REACTIONS.find(r => r.id === "init");
  const branch = M.REACTIONS.find(r => r.id === "branch");
  const coOh = M.REACTIONS.find(r => r.id === "co_oh");
  ok(init.dH > 200 && branch.dH > 50 && coOh.dH < -90,
    `entalpías clave: init=${init.dH} kJ/mol (endo), branch=${branch.dH} kJ/mol (endo), CO+OH=${coOh.dH} kJ/mol (exo)`);
}

console.log("\n── 3. Cinética de Arrhenius: metaestabilidad a 300–600 K vs autoignición a 980 K ──");
{
  const init = M.REACTIONS.find(r => r.id === "init");
  const branch = M.REACTIONS.find(r => r.id === "branch");
  const absOh = M.REACTIONS.find(r => r.id === "abs_oh");
  const pInit300 = M.arrheniusProb(init, 300);
  const pInit600 = M.arrheniusProb(init, 600);
  const pInit980 = M.arrheniusProb(init, 980);
  const pBranch300 = M.arrheniusProb(branch, 300);
  const pBranch1600 = M.arrheniusProb(branch, 1600);
  ok(pInit300 < 1e-20 && pInit600 < 1e-8 && pBranch300 < 1e-8,
    `a 300–600 K la mezcla es estable (P_init(300)=${pInit300.toExponential(1)}, P_init(600)=${pInit600.toExponential(1)})`);
  ok(pInit980 > 5e-3 && pBranch1600 > 0.1,
    `a 980 K arranca autoignición (P_init=${pInit980.toFixed(4)}) y a 1600 K la ramificación (P_branch=${pBranch1600.toFixed(3)})`);
  ok(M.arrheniusProb(absOh, 1200) > M.arrheniusProb(init, 1200),
    "abstracción por ·OH es más rápida que la iniciación térmica");
}

console.log("\n── 4. Primera Ley, estequiometría φ y llama adiabática T_ad(φ) ──");
{
  ok(near(M.equivalenceRatio(18, 36), 1.0, 1e-9), "18 CH₄ + 36 O₂ corresponde exactamente a φ = 1.00");
  const mixRich = M.mixtureCounts(1.65, 54);
  const mixLean = M.mixtureCounts(0.60, 54);
  ok(mixRich.phiActual > 1.5 && mixLean.phiActual < 0.7,
    `mixtureCounts produce mezclas coherentes (rica φ=${mixRich.phiActual.toFixed(2)}, pobre φ=${mixLean.phiActual.toFixed(2)})`);

  const TadStoich = M.adiabaticFlameTemp(18, 36, 300);
  const TadLean = M.adiabaticFlameTemp(12, 42, 300);
  const TadRich = M.adiabaticFlameTemp(26, 28, 300);
  const TadPre = M.adiabaticFlameTemp(18, 36, 980);
  ok(near(TadStoich, 2235, 30),
    `llama adiabática estequiométrica CH₄-aire: T_ad = ${Math.round(TadStoich)} K (NIST ≈ 2226 K)`);
  ok(TadStoich > TadLean && TadStoich > TadRich,
    `T_ad alcanza máximo cerca de φ = 1 (${Math.round(TadStoich)} K > pobre ${Math.round(TadLean)} K, rica ${Math.round(TadRich)} K)`);
  ok(TadPre > 2800 && TadPre < 3400,
    `mezcla precalentada a 980 K no se recorta a 2800 K (T_ad = ${Math.round(TadPre)} K)`);

  // Invariancia de trayectoria (Ley de Hess): estado intermedio conservando exactamente C=10, H=40, O=40
  const c0 = { CH4: 10, O2: 20 };
  const H0 = M.chemicalEnthalpy(c0);
  const cMid = { CH4: 4, O2: 11, CO: 6, H2O: 12 };
  const cFin = { CO2: 10, H2O: 20 };
  const T0 = M.computeTemperature(c0, H0, 0, 0, 20, 300);
  const Tmid = M.computeTemperature(cMid, H0, 0, 0, 20, 300);
  const Tfin = M.computeTemperature(cFin, H0, 0, 0, 20, 300);
  ok(near(T0, 300, 1e-9) && Tmid > T0 && Tfin > Tmid,
    `T(t) por Primera Ley crece monotónicamente al oxidar CH₄ → CO → CO₂ (${Math.round(T0)} → ${Math.round(Tmid)} → ${Math.round(Tfin)} K)`);
  ok(JSON.stringify(M.countAtoms(c0)) === JSON.stringify(M.countAtoms(cMid)) &&
     JSON.stringify(M.countAtoms(c0)) === JSON.stringify(M.countAtoms(cFin)),
    "countAtoms verifica conservación global C=10, H=40, O=40 en inicial, intermedio y final");
}

console.log("\n── 5. Regímenes de llama y geometrías moleculares 3D reales ──");
{
  ok(M.flameRegime(1.0, 300, 0, 0).label.includes("fría"), "T = 300 K sin reacciones se clasifica como mezcla fría");
  ok(!M.flameRegime(1.0, 980, 0, 0).label.includes("Llama"), "T = 980 K con rxnCount=0 no dice falsamente 'Llama'");
  ok(M.flameRegime(1.0, 600, 0, 15, 30, 6).label.includes("Ignición"), "durante ignición activa (T<750 K con radicales) dice 'Ignición en curso'");
  ok(M.flameRegime(1.0, 500, 0, 20, 0, 0).label.includes("enfriados"), "productos tras enfriarse (rxnCount>0, rad=0, T<750 K) dicen 'Productos enfriados'");
  ok(M.flameRegime(1.0, 1500, 0.65, 40, 14, 5).label.includes("azul"), "φ = 1.0 con O₂ abundante y CO transitorio sigue siendo llama azul");
  ok(M.flameRegime(1.0, 1900, 0.02, 50, 2, 1).label.includes("azul"), "φ = 1.0 caliente sin CO produce llama azul");
  ok(M.flameRegime(1.65, 1800, 0.35, 50, 0, 1).label.includes("rica"), "φ = 1.65 rica con CO produce régimen de llama rica con CO");

  const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
  const angleDeg = (a, c, b) => {
    const ux = a.x - c.x, uy = a.y - c.y, uz = a.z - c.z;
    const vx = b.x - c.x, vy = b.y - c.y, vz = b.z - c.z;
    const cos = (ux * vx + uy * vy + uz * vz) / (Math.hypot(ux, uy, uz) * Math.hypot(vx, vy, vz));
    return Math.acos(Math.max(-1, Math.min(1, cos))) * 180 / Math.PI;
  };

  const ch4 = M.MOL_3D.CH4.atoms;
  ok(near(dist(ch4[0], ch4[1]), 1.09, 1e-4) && near(angleDeg(ch4[1], ch4[0], ch4[2]), 109.47, 0.05),
    `CH₄ 3D tetraédrico exacto: d(C–H) = ${dist(ch4[0], ch4[1]).toFixed(2)} Å, ∠H-C-H = ${angleDeg(ch4[1], ch4[0], ch4[2]).toFixed(2)}°`);

  const co2 = M.MOL_3D.CO2.atoms;
  ok(near(dist(co2[0], co2[1]), 1.16, 1e-4) && near(angleDeg(co2[1], co2[0], co2[2]), 180.0, 0.01),
    `CO₂ 3D lineal exacto: d(C=O) = ${dist(co2[0], co2[1]).toFixed(2)} Å, ∠O-C-O = ${angleDeg(co2[1], co2[0], co2[2]).toFixed(1)}°`);

  const h2o = M.MOL_3D.H2O.atoms;
  ok(near(dist(h2o[0], h2o[1]), 0.96, 1e-4) && near(angleDeg(h2o[1], h2o[0], h2o[2]), 104.5, 0.05),
    `H₂O 3D angular exacto: d(O–H) = ${dist(h2o[0], h2o[1]).toFixed(2)} Å, ∠H-O-H = ${angleDeg(h2o[1], h2o[0], h2o[2]).toFixed(1)}°`);
}

console.log(fail ? `\n${fail} prueba(s) fallaron\n` : "\ntodas las pruebas pasaron\n");
process.exit(fail ? 1 : 0);
