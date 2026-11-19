/* combustion.js — Simulador de combustión del metano: CH₄ + 2 O₂ → CO₂ + 2 H₂O.
 * Mecanismo reducido didáctico en cadena ramificada (29 reacciones bimoleculares +
 * disociación térmica de ·CHO = 30 pasos irreversibles) con:
 *   - Conservación exacta de átomos (C, H, O) y barrera Eₐ ≥ ΔH° en cada paso endotérmico
 *   - Balance energético de Primera Ley con ΔH_f°(298 K) NIST/ATcT y Cp medio en llama
 *     (aire seco 3.76 N₂/O₂ → T_ad aprox. ≈ 2235 K sin disociación de equilibrio)
 *   - Relación de equivalencia φ (pobre, estequiométrica azul y rica con CO tóxico + H₂)
 *   - Chispa eléctrica local (núcleo caliente de radicales ·O· sin precalentar toda la cámara)
 *   - Visor 3D interactivo de geometrías moleculares (ratón, táctil y teclado)
 * Funciones puras exportadas para verificación en Node (tests/check.js). */
(function () {
"use strict";

// ==================== CONSTANTES Y TERMOQUÍMICA ====================
const R_GAS = 8.3144626e-3;   // kJ/(mol·K)
const N2_PER_O2 = 3.76;       // relación molar N₂/O₂ del aire seco
const CP_N2 = 35.5;           // J/(mol·K) medio de N₂ a temperaturas de llama (~1400 K)
const DH_VAP_H2O = 44.0;      // kJ/mol de vaporización del agua a 298 K

// Especies químicas: (C, H, O), ΔH_f°(298 K) [kJ/mol], Cp medio en llama [J/(mol·K)], masa [g/mol]
const SPECIES = {
  CH4:  { label:"CH₄",  name:"Metano",                C:1, H:4, O:0, Hf:-74.6,  Cp:72, M:16.04, rad:false, r:9.5, col:"#6e7681" },
  O2:   { label:"O₂",   name:"Oxígeno molecular",     C:0, H:0, O:2, Hf:0.0,    Cp:36, M:32.00, rad:false, r:8.2, col:"#f85149" },
  H2:   { label:"H₂",   name:"Hidrógeno molecular",   C:0, H:2, O:0, Hf:0.0,    Cp:31, M:2.02,  rad:false, r:5.5, col:"#e6edf3" },
  CH3:  { label:"·CH₃", name:"Radical metilo",        C:1, H:3, O:0, Hf:146.7,  Cp:60, M:15.03, rad:true,  r:8.4, col:"#b197fc" },
  OH:   { label:"·OH",  name:"Radical hidroxilo",     C:0, H:1, O:1, Hf:37.3,   Cp:33, M:17.01, rad:true,  r:6.8, col:"#ff8787" },
  O:    { label:"·O·",  name:"Oxígeno atómico",       C:0, H:0, O:1, Hf:249.2,  Cp:21, M:16.00, rad:true,  r:5.8, col:"#ff6b6b" },
  H:    { label:"·H",   name:"Hidrógeno atómico",     C:0, H:1, O:0, Hf:218.0,  Cp:21, M:1.01,  rad:true,  r:4.4, col:"#ffd43b" },
  HO2:  { label:"HO₂·", name:"Radical hidroperoxilo", C:0, H:1, O:2, Hf:12.3,   Cp:48, M:33.01, rad:true,  r:8.2, col:"#da77f2" },
  CH2O: { label:"CH₂O", name:"Formaldehído",          C:1, H:2, O:1, Hf:-109.2, Cp:58, M:30.03, rad:false, r:8.6, col:"#c084fc" },
  CHO:  { label:"·CHO", name:"Radical formilo",       C:1, H:1, O:1, Hf:42.5,   Cp:48, M:29.02, rad:true,  r:8.0, col:"#f783ac" },
  CO:   { label:"CO",   name:"Monóxido de carbono",   C:1, H:0, O:1, Hf:-110.5, Cp:35, M:28.01, rad:false, r:7.8, col:"#ffa657" },
  CO2:  { label:"CO₂",  name:"Dióxido de carbono",    C:1, H:0, O:2, Hf:-393.5, Cp:56, M:44.01, rad:false, r:9.4, col:"#3fb950" },
  H2O:  { label:"H₂O",  name:"Vapor de agua",         C:0, H:2, O:1, Hf:-241.8, Cp:46, M:18.02, rad:false, r:7.6, col:"#58a6ff" }
};

// Energías medias de disociación de enlace (kJ/mol, Atkins/NIST; C=O específico de CO₂; primer BDE CH₃–H ≈ 439)
const BOND_ENERGIES = { CH: 413, CH_FIRST: 439, OO: 498, CO_CO2: 799, OH: 463 };

function bondBalance() {
  const brokenCH = 4 * BOND_ENERGIES.CH;
  const brokenOO = 2 * BOND_ENERGIES.OO;
  const broken = brokenCH + brokenOO;
  const formedCO = 2 * BOND_ENERGIES.CO_CO2;
  const formedOH = 4 * BOND_ENERGIES.OH;
  const formed = formedCO + formedOH;
  return { brokenCH, brokenOO, broken, formedCO, formedOH, formed, net: broken - formed };
}

function hessReactionEnthalpy(liquidWater = false) {
  const gas = SPECIES.CO2.Hf + 2 * SPECIES.H2O.Hf - SPECIES.CH4.Hf - 2 * SPECIES.O2.Hf;
  return liquidWater ? gas - 2 * DH_VAP_H2O : gas;
}

// ==================== MECANISMO REDUCIDO EN CADENA RAMIFICADA ====================
// Todas las reacciones conservan exactamente (C, H, O) y cumplen Eₐ ≥ ΔH° en pasos endotérmicos.
// P0 es el factor preexponencial efectivo por colisión escalado a la frecuencia de choque de la cámara 2D
// (p. ej. init mantiene Eₐ real = 238 kJ/mol con P0 = 3.5e10 y techo estérico maxP = 0.04 → nulo a 300–600 K, autoignición real a ~980 K).
const REACTIONS = [
  { id:"init",    r:["CH4","O2"],  p:["CH3","HO2"],    Ea:238, P0:3.5e10, maxP:0.04, cat:"iniciación",    eq:"CH₄ + O₂ → ·CH₃ + HO₂·" },
  { id:"abs_oh",  r:["CH4","OH"],  p:["CH3","H2O"],    Ea:14,  P0:4.0,    cat:"propagación",   eq:"CH₄ + ·OH → ·CH₃ + H₂O" },
  { id:"abs_o",   r:["CH4","O"],   p:["CH3","OH"],     Ea:36,  P0:22.0,   cat:"propagación",   eq:"CH₄ + ·O· → ·CH₃ + ·OH" },
  { id:"abs_h",   r:["CH4","H"],   p:["CH3","H2"],     Ea:45,  P0:15.0,   cat:"propagación",   eq:"CH₄ + ·H → ·CH₃ + H₂" },
  { id:"branch",  r:["H","O2"],    p:["OH","O"],       Ea:70,  P0:600.0,  cat:"ramificación",  eq:"·H + O₂ → ·OH + ·O·" },
  { id:"h2_oh",   r:["H2","OH"],   p:["H2O","H"],      Ea:18,  P0:4.0,    cat:"propagación",   eq:"H₂ + ·OH → H₂O + ·H" },
  { id:"h2_o",    r:["H2","O"],    p:["OH","H"],       Ea:38,  P0:16.0,   cat:"ramificación",  eq:"H₂ + ·O· → ·OH + ·H" },
  { id:"ch3_o",   r:["CH3","O"],   p:["CH2O","H"],     Ea:4,   P0:2.0,    cat:"oxidación",     eq:"·CH₃ + ·O· → CH₂O + ·H" },
  { id:"ch3_o2",  r:["CH3","O2"],  p:["CH2O","OH"],    Ea:38,  P0:38.0,   cat:"oxidación",     eq:"·CH₃ + O₂ → CH₂O + ·OH" },
  // Paso efectivo condensado de propagación (GRI-Mech: ·CH₃ + HO₂· → CH₃O· + ·OH → CH₂O + ·OH + ·H, ΔH = −12.9 kJ/mol)
  { id:"ch3_ho2", r:["CH3","HO2"], p:["CH2O","OH","H"], Ea:8,  P0:2.0,    cat:"propagación",   eq:"·CH₃ + HO₂· → CH₂O + ·OH + ·H (vía CH₃O·)" },
  { id:"ch2o_oh", r:["CH2O","OH"], p:["CHO","H2O"],    Ea:5,   P0:2.5,    cat:"oxidación",     eq:"CH₂O + ·OH → ·CHO + H₂O" },
  { id:"ch2o_o",  r:["CH2O","O"],  p:["CHO","OH"],     Ea:15,  P0:4.0,    cat:"oxidación",     eq:"CH₂O + ·O· → ·CHO + ·OH" },
  { id:"ch2o_h",  r:["CH2O","H"],  p:["CHO","H2"],     Ea:16,  P0:3.5,    cat:"oxidación",     eq:"CH₂O + ·H → ·CHO + H₂" },
  { id:"cho_o2",  r:["CHO","O2"],  p:["CO","HO2"],     Ea:4,   P0:1.5,    cat:"oxidación",     eq:"·CHO + O₂ → CO + HO₂·" },
  { id:"cho_oh",  r:["CHO","OH"],  p:["CO","H2O"],     Ea:2,   P0:1.0,    cat:"oxidación",     eq:"·CHO + ·OH → CO + H₂O" },
  { id:"cho_o",   r:["CHO","O"],   p:["CO","OH"],      Ea:2,   P0:1.2,    cat:"oxidación",     eq:"·CHO + ·O· → CO + ·OH" },
  { id:"cho_h",   r:["CHO","H"],   p:["CO","H2"],      Ea:2,   P0:1.0,    cat:"oxidación",     eq:"·CHO + ·H → CO + H₂" },
  { id:"co_oh",   r:["CO","OH"],   p:["CO2","H"],      Ea:12,  P0:6.0,    cat:"cierre CO₂",    eq:"CO + ·OH → CO₂ + ·H" },
  { id:"co_ho2",  r:["CO","HO2"],  p:["CO2","OH"],     Ea:62,  P0:320.0,  cat:"cierre CO₂",    eq:"CO + HO₂· → CO₂ + ·OH" },
  { id:"co_o",    r:["CO","O"],    p:["CO2"],          Ea:18,  P0:2.2,    cat:"cierre CO₂",    eq:"CO + ·O· (+M) → CO₂" },
  { id:"ho2_h",   r:["HO2","H"],   p:["OH","OH"],      Ea:4,   P0:1.5,    cat:"ramificación",  eq:"HO₂· + ·H → 2 ·OH" },
  { id:"ho2_o",   r:["HO2","O"],   p:["OH","O2"],      Ea:2,   P0:1.2,    cat:"propagación",   eq:"HO₂· + ·O· → ·OH + O₂" },
  // Paso efectivo condensado (2 HO₂· → H₂O₂ + O₂ → 2 ·OH + O₂ en alta T; Eₐ = 54 ≥ ΔH = +50.0)
  { id:"ho2_ho2", r:["HO2","HO2"], p:["OH","OH","O2"], Ea:54,  P0:220.0,  cat:"efectivo",      eq:"2 HO₂· → 2 ·OH + O₂ (vía H₂O₂)" },
  { id:"ho2_oh",  r:["HO2","OH"],  p:["H2O","O2"],     Ea:2,   P0:0.16,   cat:"terminación",   eq:"HO₂· + ·OH → H₂O + O₂" },
  { id:"oh_oh",   r:["OH","OH"],   p:["H2O","O"],      Ea:8,   P0:0.75,   cat:"propagación",   eq:"2 ·OH → H₂O + ·O·" },
  // Recombinaciones trimoleculares (+M): cierran el pool de radicales liberando su entalpía de formación
  { id:"h_oh",    r:["H","OH"],    p:["H2O"],          Ea:0,   P0:0.14,   cat:"terminación",   eq:"·H + ·OH (+M) → H₂O" },
  { id:"h_h",     r:["H","H"],     p:["H2"],           Ea:0,   P0:0.10,   cat:"terminación",   eq:"2 ·H (+M) → H₂" },
  { id:"o_o",     r:["O","O"],     p:["O2"],           Ea:0,   P0:0.08,   cat:"terminación",   eq:"2 ·O· (+M) → O₂" },
  { id:"ch3_h",   r:["CH3","H"],   p:["CH4"],          Ea:0,   P0:0.04,   cat:"terminación",   eq:"·CH₃ + ·H (+M) → CH₄" }
];

const UNIMOL_CHO = {
  id:"cho_m", r:["CHO"], p:["CO","H"], Ea:68, P0:800.0, cat:"disociación", eq:"·CHO (+M) → CO + ·H",
  dH: +(SPECIES.CO.Hf + SPECIES.H.Hf - SPECIES.CHO.Hf).toFixed(2)
};

// Cada par no ordenado de especies tiene un único canal representativo en este mecanismo reducido
const RXN_MAP = new Map();
for (const rx of REACTIONS) {
  const sum = (arr, k) => arr.reduce((acc, s) => acc + SPECIES[s][k], 0);
  rx.dH = +(sum(rx.p, "Hf") - sum(rx.r, "Hf")).toFixed(2);
  const key = rx.r.slice().sort().join("+");
  RXN_MAP.set(key, rx);
}

/** Probabilidad de Arrhenius efectiva por colisión: P = min(maxP, P0 · exp(−Eₐ / (R·T))) */
function arrheniusProb(rx, T) {
  if (!(T > 1)) return 0;
  return Math.min(rx.maxP || 1, rx.P0 * Math.exp(-rx.Ea / (R_GAS * T)));
}

/** Relación de equivalencia φ = (n_CH4 / n_O2) / (1/2) = 2·n_CH4 / n_O2 */
function equivalenceRatio(nCH4, nO2) {
  if (!(nO2 > 0)) return Infinity;
  return (2 * nCH4) / nO2;
}

/** Reparte N moléculas iniciales en (nCH4, nO2) según φ ∈ [0.5, 2.0] */
function mixtureCounts(phi, totalMols = 54) {
  const f = Math.max(0.4, Math.min(2.2, phi));
  const nCH4 = Math.max(4, Math.round(totalMols * f / (f + 2)));
  const nO2 = Math.max(6, totalMols - nCH4);
  return { nCH4, nO2, phiActual: equivalenceRatio(nCH4, nO2) };
}

/** Conteo de átomos totales (C, H, O) de un diccionario de conteos por especie */
function countAtoms(counts) {
  let C = 0, H = 0, O = 0;
  for (const [k, n] of Object.entries(counts)) {
    const sp = SPECIES[k];
    if (!sp || !n) continue;
    C += sp.C * n; H += sp.H * n; O += sp.O * n;
  }
  return { C, H, O };
}

/** Entalpía química total Σ n_i ΔH_f°(i) en kJ (escala 1 partícula = 1 mol) */
function chemicalEnthalpy(counts) {
  let H = 0;
  for (const [k, n] of Object.entries(counts)) {
    if (SPECIES[k] && n) H += SPECIES[k].Hf * n;
  }
  return H;
}

/** Capacidad calorífica isobárica media Cp [kJ/K] de la mezcla + N₂ del aire (3.76·nO2_0) */
function totalHeatCapacity(counts, nO2Init) {
  let cpJ = nO2Init * N2_PER_O2 * CP_N2;
  for (const [k, n] of Object.entries(counts)) {
    if (SPECIES[k] && n) cpJ += SPECIES[k].Cp * n;
  }
  return cpJ * 1e-3; // kJ/K
}

/** Temperatura por balance de entalpía (Cp medio): T = T₀ + (H_chem(0) − H_chem(t) + Q_spark − Q_loss) / Cp */
function computeTemperature(counts, Hchem0, sparkKJ, lossKJ, nO2Init, T0 = 300) {
  const Hnow = chemicalEnthalpy(counts);
  const dQ = (Hchem0 - Hnow) + sparkKJ - lossKJ;
  const Cp = Math.max(0.1, totalHeatCapacity(counts, nO2Init));
  return Math.max(50, Math.min(4000, T0 + dQ / Cp));
}

/** Temperatura límite adiabática aproximada T_ad(φ) en aire seco (Cp medio, sin disociación de equilibrio).
 *  Para mezcla rica (φ > 1) particiona los productos entre CO, CO₂, H₂ y H₂O conservando exactamente C, H y O. */
function adiabaticFlameTemp(nCH4, nO2, T0 = 300) {
  const init = { CH4: nCH4, O2: nO2 };
  const H0 = chemicalEnthalpy(init);
  let fin;
  if (nO2 >= 2 * nCH4) {
    fin = { CO2: nCH4, H2O: 2 * nCH4, O2: nO2 - 2 * nCH4 };
  } else {
    // Oxidación parcial (CO + H₂) + oxidación adicional del O₂ sobrante hacia H₂O y CO₂ (conservando C, H, O)
    const burned = Math.min(nCH4, 2 * nO2); // CH₄ + 0.5 O₂ → CO + 2 H₂
    let oRem = 2 * nO2 - burned;            // átomos de O restantes tras formar `burned` CO
    const nH2O = Math.min(2 * burned, oRem * 0.65);
    oRem -= nH2O;
    const nCO2 = Math.min(burned, oRem);
    oRem -= nCO2;
    const extraH2O = Math.min(2 * burned - nH2O, oRem);
    const totH2O = nH2O + extraH2O;
    fin = {
      CH4: nCH4 - burned,
      CO: burned - nCO2,
      CO2: nCO2,
      H2O: totH2O,
      H2: 2 * burned - totH2O
    };
  }
  return computeTemperature(fin, H0, 0, 0, nO2, T0);
}

/** Estado y color cualitativo de la llama según régimen y avance real de reacción */
function flameRegime(phi, T, coRatio, rxnCount = 1, o2Left = 0, radCount = 0) {
  if (rxnCount === 0) {
    if (T >= 850) return { label:"Mezcla caliente · retraso de autoignición…", col:"#f0c000", glow:"rgba(240,192,0,0.08)" };
    return { label:"Mezcla fría (metaestable)", col:"#58a6ff", glow:"rgba(88,166,255,0.04)" };
  }
  if (T < 750) {
    if (radCount > 0) return { label:"Ignición en curso · cascada de radicales…", col:"#38bdf8", glow:"rgba(56,189,248,0.10)" };
    return { label:"Productos enfriados", col:"#79c0ff", glow:"rgba(121,192,255,0.04)" };
  }
  if (phi > 1.25 || (o2Left <= 3 && coRatio > 0.25)) {
    return { label:"Llama rica · CO tóxico (hollín no simulado)", col:"#ffa657", glow:"rgba(255,166,87,0.16)" };
  }
  if (phi < 0.85) return { label:"Llama azul pálido · mezcla pobre", col:"#79c0ff", glow:"rgba(88,166,255,0.13)" };
  return { label:"Llama azul limpia · combustión completa", col:"#38bdf8", glow:"rgba(56,189,248,0.16)" };
}

// ==================== GEOMETRÍAS 3D REALES (Å) ====================
const TET = 1 / Math.sqrt(3);
const H2O_HALF = (104.5 / 2) * Math.PI / 180;
const CH2O_HALF = (116.5 / 2) * Math.PI / 180;

const MOL_3D = {
  CH4: {
    formula: "CH₄", name: "Metano", geom: "Tetraédrica · sp³ · 109.5°", bondInfo: "4 C–H (109 pm)",
    atoms: [
      { el:"C", x:0, y:0, z:0 },
      { el:"H", x: 1.09*TET, y: 1.09*TET, z: 1.09*TET },
      { el:"H", x: 1.09*TET, y:-1.09*TET, z:-1.09*TET },
      { el:"H", x:-1.09*TET, y: 1.09*TET, z:-1.09*TET },
      { el:"H", x:-1.09*TET, y:-1.09*TET, z: 1.09*TET }
    ],
    bonds: [[0,1,1],[0,2,1],[0,3,1],[0,4,1]]
  },
  O2: {
    formula: "O₂", name: "Oxígeno", geom: "Lineal · doble enlace", bondInfo: "O=O (121 pm)",
    atoms: [
      { el:"O", x:-0.605, y:0, z:0 },
      { el:"O", x: 0.605, y:0, z:0 }
    ],
    bonds: [[0,1,2]]
  },
  CH2O: {
    formula: "CH₂O", name: "Formaldehído", geom: "Trigonal plana · sp² · 116.5°", bondInfo: "C=O 121 · C–H 111 pm",
    atoms: [
      { el:"C", x:-0.15, y:0, z:0 },
      { el:"O", x: 1.06, y:0, z:0 },
      { el:"H", x:-0.15 - 1.11*Math.cos(CH2O_HALF), y: 1.11*Math.sin(CH2O_HALF), z:0 },
      { el:"H", x:-0.15 - 1.11*Math.cos(CH2O_HALF), y:-1.11*Math.sin(CH2O_HALF), z:0 }
    ],
    bonds: [[0,1,2],[0,2,1],[0,3,1]]
  },
  CO: {
    formula: "CO", name: "Monóxido de carbono", geom: "Lineal · triple enlace", bondInfo: "C≡O (113 pm · tóxico)",
    atoms: [
      { el:"C", x:-0.565, y:0, z:0 },
      { el:"O", x: 0.565, y:0, z:0 }
    ],
    bonds: [[0,1,3]]
  },
  CO2: {
    formula: "CO₂", name: "Dióxido de carbono", geom: "Lineal · sp · 180.0°", bondInfo: "2 C=O (116 pm)",
    atoms: [
      { el:"C", x:0, y:0, z:0 },
      { el:"O", x:-1.16, y:0, z:0 },
      { el:"O", x: 1.16, y:0, z:0 }
    ],
    bonds: [[0,1,2],[0,2,2]]
  },
  H2O: {
    formula: "H₂O", name: "Agua (gas)", geom: "Angular · 104.5°", bondInfo: "2 O–H (96 pm)",
    atoms: [
      { el:"O", x:0, y: 0.25, z:0 },
      { el:"H", x:-0.96*Math.sin(H2O_HALF), y:0.25 - 0.96*Math.cos(H2O_HALF), z:0 },
      { el:"H", x: 0.96*Math.sin(H2O_HALF), y:0.25 - 0.96*Math.cos(H2O_HALF), z:0 }
    ],
    bonds: [[0,1,1],[0,2,1]]
  },
  OH: {
    formula: "·OH", name: "Radical hidroxilo", geom: "1 e⁻ desapareado", bondInfo: "O–H (97 pm · propagador)",
    atoms: [
      { el:"O", x:-0.32, y:0, z:0, radical:true },
      { el:"H", x: 0.65, y:0, z:0 }
    ],
    bonds: [[0,1,1]]
  }
};

// ==================== PRESETS DE MEZCLA ====================
const PRESETS = [
  { id:"stoich", label:"Estequiométrica", sub:"φ = 1.0 (azul)", phi:1.0, T0:300,
    desc:"Proporción exacta 1 CH₄ : 2 O₂. Fría a 300 K es metaestable; una chispa local crea radicales ·O· que encienden la cadena hacia CO₂ y H₂O (T_ad ≈ 2235 K)." },
  { id:"rich",   label:"Rica (poco O₂)", sub:"φ = 1.65 (CO)", phi:1.65, T0:300,
    desc:"Falta oxígeno (φ > 1.25): la combustión incompleta acumula monóxido de carbono (CO) tóxico e hidrógeno (H₂). En llamas reales el exceso de carbono genera hollín amarillo." },
  { id:"lean",   label:"Pobre (exceso O₂)", sub:"φ = 0.60", phi:0.60, T0:300,
    desc:"Exceso de aire: todo el metano pasa a CO₂ + H₂O con temperatura máxima menor porque el O₂ y N₂ sobrantes absorben parte del calor." },
  { id:"auto",   label:"Autoignición", sub:"T₀ = 980 K", phi:1.0, T0:980,
    desc:"A 980 K la cola térmica supera la barrera de iniciación (Eₐ = 238 kJ/mol) sin chispa externa: tras un breve retraso de ignición estalla la cadena ramificada." }
];

// ==================== INTERFAZ DE NAVEGADOR ====================
if (typeof document !== "undefined") (function () {
const $ = id => document.getElementById(id);

const cvReactor = $("reactor"), ctxR = cvReactor.getContext("2d");
const cvKin = $("kinCanvas"), ctxK = cvKin.getContext("2d");
const cvProf = $("profCanvas"), ctxP = cvProf.getContext("2d");
const cvMol = $("molCanvas"), ctxM = cvMol.getContext("2d");

const TAU = Math.PI * 2;
const reducedMotion = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

function roundRectPath(ctx, x, y, w, h, r) {
  if (typeof ctx.roundRect === "function") { ctx.roundRect(x, y, w, h, r); return; }
  const rr = Math.min(r, w / 2, h / 2);
  ctx.moveTo(x + rr, y);
  ctx.arcTo(x + w, y, x + w, y + h, rr);
  ctx.arcTo(x + w, y + h, x, y + h, rr);
  ctx.arcTo(x, y + h, x, y, rr);
  ctx.arcTo(x, y, x + w, y, rr);
  ctx.closePath();
}

function fitCanvas(cv, ctx) {
  const dpr = window.devicePixelRatio || 1;
  const W = cv.clientWidth, H = cv.clientHeight;
  if (!W || !H) return { W: 0, H: 0 };
  const rw = Math.round(W * dpr), rh = Math.round(H * dpr);
  if (cv.width !== rw || cv.height !== rh) { cv.width = rw; cv.height = rh; }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { W, H };
}

const state = {
  presetId: "stoich",
  phi: 1.0,
  T0: 300,
  T: 300,
  adiabatic: true,
  paused: false,
  speed: reducedMotion ? 0.5 : 1,
  simTime: 0,
  parts: [],
  flashes: [],
  sparkFX: null,
  history: [],
  initCH4: 18,
  initO2: 36,
  initAtoms: { C: 18, H: 72, O: 72 },
  Hchem0: 0,
  sparkKJ: 0,
  sparkCount: 0,
  lossKJ: 0,
  rxnCount: 0,
  molKey: "CH4",
  rotY: 0.45,
  rotX: 0.30,
  molDirty: true,
  reactorDirty: true,
  infoH: 18,
  W: 0, H: 0
};

let pid = 0;
function currentCounts() {
  const c = {};
  for (const p of state.parts) c[p.sp] = (c[p.sp] || 0) + 1;
  return c;
}

function makeMol(sp, x, y, vx, vy, hot = 0) {
  const s = SPECIES[sp];
  return {
    id: ++pid, sp, x, y, vx, vy, hot,
    r: s.r, M: s.M, ang: Math.random() * TAU, spin: (Math.random() - 0.5) * 2.6
  };
}

function thermalSpeed(sp, T) {
  const base = 155 * Math.sqrt((Math.max(T, 200) / 300) * (16 / SPECIES[sp].M));
  return Math.min(560, base * (0.68 + Math.random() * 0.64));
}

function rescaleVelocities(Told, Tnew) {
  const f = Math.sqrt(Math.max(Tnew, 200) / Math.max(Told, 200));
  for (const p of state.parts) {
    p.vx *= f; p.vy *= f;
    const sp = Math.hypot(p.vx, p.vy);
    const vmax = Math.min(560, 240 * Math.sqrt(Math.max(Tnew, 200) / 300) * Math.sqrt(16 / p.M));
    if (sp > vmax) { p.vx *= vmax / sp; p.vy *= vmax / sp; }
  }
}

function resetReactor() {
  const sz = fitCanvas(cvReactor, ctxR);
  state.W = Math.max(sz.W, 240);
  state.H = Math.max(sz.H, 180);
  const mix = mixtureCounts(state.phi, 54);
  state.initCH4 = mix.nCH4;
  state.initO2 = mix.nO2;
  state.parts = [];
  state.flashes = [];
  state.sparkFX = null;
  state.simTime = 0;
  state.T = state.T0;
  state.sparkKJ = 0;
  state.sparkCount = 0;
  state.lossKJ = 0;
  state.rxnCount = 0;
  pid = 0;
  histAccum = 0;
  uiAccum = 0;
  $("btnSpark").disabled = false;

  for (let i = 0; i < mix.nCH4; i++) {
    const a = Math.random() * TAU, v = thermalSpeed("CH4", state.T0);
    state.parts.push(makeMol("CH4", 24 + Math.random() * (state.W - 48), 24 + Math.random() * (state.H - 48), Math.cos(a) * v, Math.sin(a) * v));
  }
  for (let i = 0; i < mix.nO2; i++) {
    const a = Math.random() * TAU, v = thermalSpeed("O2", state.T0);
    state.parts.push(makeMol("O2", 24 + Math.random() * (state.W - 48), 24 + Math.random() * (state.H - 48), Math.cos(a) * v, Math.sin(a) * v));
  }

  const c0 = currentCounts();
  state.initAtoms = countAtoms(c0);
  state.Hchem0 = chemicalEnthalpy(c0);
  state.history = [];
  recordHistory(c0);
  clearRxnLog();
  state.reactorDirty = true;
  state.molDirty = true;
  updateSideUI(c0);
  drawKinetics();
  drawProfile(c0);
  drawMol3D();
}

function recordHistory(c) {
  const nRad = (c.CH3 || 0) + (c.OH || 0) + (c.O || 0) + (c.H || 0) + (c.HO2 || 0) + (c.CHO || 0);
  state.history.push({
    t: state.simTime,
    T: state.T,
    CH4: c.CH4 || 0,
    O2: c.O2 || 0,
    CO2: c.CO2 || 0,
    H2O: c.H2O || 0,
    CO: c.CO || 0,
    CH2O: c.CH2O || 0,
    rad: nRad
  });
  if (state.history.length > 260) state.history.shift();
}

function clearRxnLog() {
  const el = $("rxnLog");
  el.textContent = "";
  const n = document.createElement("div");
  n.className = "note";
  n.textContent = state.T0 >= 850
    ? `A T₀ = ${state.T0} K la mezcla supera el umbral térmico: observa el retraso de autoignición sin chispa.`
    : `A T₀ = ${state.T0} K el metano y el oxígeno chocan sin reaccionar por la barrera de activación (Eₐ = 238 kJ/mol). Pulsa ⚡ Chispa de ignición.`;
  el.appendChild(n);
}

function appendRxnLog(rx) {
  const el = $("rxnLog");
  const note = el.querySelector(".note");
  if (note) el.textContent = "";
  const row = document.createElement("div");
  row.className = "rlog-row";
  const eq = document.createElement("span");
  eq.className = "rlog-eq";
  eq.textContent = rx.eq;
  const meta = document.createElement("span");
  meta.className = rx.dH < 0 ? "rlog-exo" : "rlog-endo";
  meta.textContent = `${rx.dH > 0 ? "+" : ""}${rx.dH.toFixed(0)} kJ/mol`;
  row.append(eq, meta);
  el.prepend(row);
  while (el.children.length > 35) el.lastChild.remove();
}

// ---------- Chispa eléctrica local (núcleo de radicales ·O· sin precalentar toda la cámara) ----------
function fireSpark() {
  setPaused(false);
  const cx = state.W * 0.5, cy = state.H * 0.5;
  const o2List = state.parts
    .filter(p => p.sp === "O2")
    .sort((a, b) => Math.hypot(a.x - cx, a.y - cy) - Math.hypot(b.x - cx, b.y - cy));
  const toSplit = o2List.slice(0, 3);
  if (!toSplit.length) {
    appendRxnLog({ eq: "⚡ Sin O₂ libre para disociar en la chispa", dH: 0 });
    return;
  }
  state.sparkFX = { x: cx, y: cy, age: 0 };
  const kill = new Set(toSplit);
  const added = [];
  for (const m of toSplit) {
    const a = Math.random() * TAU, kick = 220;
    added.push(makeMol("O", m.x + 4, m.y, m.vx + Math.cos(a) * kick, m.vy + Math.sin(a) * kick, 1850));
    added.push(makeMol("O", m.x - 4, m.y, m.vx - Math.cos(a) * kick, m.vy - Math.sin(a) * kick, 1850));
  }
  state.parts = state.parts.filter(p => !kill.has(p));
  for (const a of added) state.parts.push(a);

  // Calienta únicamente el núcleo local cercano a la chispa (plasma local ~1950 K)
  for (const p of state.parts) {
    if (Math.hypot(p.x - cx, p.y - cy) < 200) p.hot = Math.max(p.hot, 1950);
  }
  const cNow = currentCounts();
  const Cp = totalHeatCapacity(cNow, state.initO2);
  const dissocKJ = toSplit.length * (2 * SPECIES.O.Hf);
  // Energía acotada de chispa: disociación O=O → 2 ·O· + mínimo plasma térmico (+80 K en 1ª chispa)
  const localHeat = state.sparkCount === 0 ? 80 * Cp : 0;
  state.sparkKJ += dissocKJ + localHeat;
  state.sparkCount++;
  const Told = state.T;
  state.T = computeTemperature(cNow, state.Hchem0, state.sparkKJ, state.lossKJ, state.initO2, state.T0);
  rescaleVelocities(Told, state.T);
  if (state.sparkCount >= 2) $("btnSpark").disabled = true;
  appendRxnLog({ eq: `⚡ Chispa local: ${toSplit.length} O₂ → ${toSplit.length * 2} ·O·`, dH: dissocKJ });
  state.reactorDirty = true;
  updateSideUI(cNow);
  drawKinetics();
  drawProfile(cNow);
}

// ---------- Dinámica 2D (2 subpasos/frame) y química en cadena ----------
let histAccum = 0, uiAccum = 0;
function stepSub(dt) {
  state.simTime += dt;
  const W = state.W, H = state.H;

  if (!state.adiabatic && state.T > state.T0 + 1) {
    const cNow = currentCounts();
    const Cp = totalHeatCapacity(cNow, state.initO2);
    state.lossKJ += 0.09 * (state.T - state.T0) * Cp * dt;
    const Told = state.T;
    state.T = computeTemperature(cNow, state.Hchem0, state.sparkKJ, state.lossKJ, state.initO2, state.T0);
    rescaleVelocities(Told, state.T);
  }

  for (const p of state.parts) {
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    p.ang += p.spin * dt;
    if (p.hot > state.T) p.hot = Math.max(state.T, p.hot - 180 * dt);
    if (p.x < p.r) { p.x = p.r; p.vx = Math.abs(p.vx); }
    else if (p.x > W - p.r) { p.x = W - p.r; p.vx = -Math.abs(p.vx); }
    if (p.y < p.r) { p.y = p.r; p.vy = Math.abs(p.vy); }
    else if (p.y > H - p.r) { p.y = H - p.r; p.vy = -Math.abs(p.vy); }
  }

  const dead = new Set();
  const spawned = [];
  let reacted = false;

  // Disociación térmica unimolecular de ·CHO (+M) → CO + ·H
  for (const p of state.parts) {
    if (p.sp === "CHO") {
      const teff = Math.max(state.T, p.hot);
      if (Math.random() < arrheniusProb(UNIMOL_CHO, teff) * Math.min(1, dt * 6)) {
        dead.add(p);
        const a = Math.random() * TAU;
        spawned.push(makeMol("CO", p.x + 3, p.y, p.vx + Math.cos(a) * 35, p.vy + Math.sin(a) * 35, teff));
        spawned.push(makeMol("H", p.x - 3, p.y, p.vx - Math.cos(a) * 110, p.vy - Math.sin(a) * 110, teff));
        state.flashes.push({ x: p.x, y: p.y, age: 0, col: "#ffa657", r: 16 });
        state.rxnCount++;
        reacted = true;
        appendRxnLog(UNIMOL_CHO);
      }
    }
  }

  // Colisiones bimoleculares elásticas + reacciones químicas
  const parts = state.parts;
  for (let i = 0; i < parts.length; i++) {
    const a = parts[i];
    if (dead.has(a)) continue;
    for (let j = i + 1; j < parts.length; j++) {
      const b = parts[j];
      if (dead.has(b)) continue;
      const dx = b.x - a.x, dy = b.y - a.y;
      const rReact = a.r + b.r + 12.5;
      const d2 = dx * dx + dy * dy;
      if (d2 >= rReact * rReact) continue;

      const key = a.sp < b.sp ? `${a.sp}+${b.sp}` : `${b.sp}+${a.sp}`;
      const rx = RXN_MAP.get(key);
      const teff = Math.max(state.T, a.hot, b.hot);
      if (rx && Math.random() < arrheniusProb(rx, teff)) {
        dead.add(a); dead.add(b);
        const mx = (a.x + b.x) * 0.5, my = (a.y + b.y) * 0.5;
        const Mtot = a.M + b.M;
        const vcx = (a.M * a.vx + b.M * b.vx) / Mtot;
        const vcy = (a.M * a.vy + b.M * b.vy) / Mtot;
        const baseAng = Math.random() * TAU;
        const kickBase = rx.dH < 0 ? 58 : 22;
        const prodHot = rx.dH < -40 ? Math.max(teff, 1650) : teff;

        if (rx.p.length === 1) {
          spawned.push(makeMol(rx.p[0], mx, my, vcx, vcy, prodHot));
        } else if (rx.p.length === 2) {
          const m1 = SPECIES[rx.p[0]].M, m2 = SPECIES[rx.p[1]].M;
          const kx = Math.cos(baseAng) * kickBase, ky = Math.sin(baseAng) * kickBase;
          spawned.push(makeMol(rx.p[0], mx + 3, my + 2, vcx + kx * (m2 / (m1 + m2)), vcy + ky * (m2 / (m1 + m2)), prodHot));
          spawned.push(makeMol(rx.p[1], mx - 3, my - 2, vcx - kx * (m1 / (m1 + m2)), vcy - ky * (m1 / (m1 + m2)), prodHot));
        } else {
          const m0 = SPECIES[rx.p[0]].M, m1 = SPECIES[rx.p[1]].M, m2 = SPECIES[rx.p[2]].M;
          const v0x = Math.cos(baseAng) * kickBase, v0y = Math.sin(baseAng) * kickBase;
          const v1x = Math.cos(baseAng + 2.1) * kickBase, v1y = Math.sin(baseAng + 2.1) * kickBase;
          const v2x = -(m0 * v0x + m1 * v1x) / m2, v2y = -(m0 * v0y + m1 * v1y) / m2;
          spawned.push(makeMol(rx.p[0], mx + 3, my, vcx + v0x, vcy + v0y, prodHot));
          spawned.push(makeMol(rx.p[1], mx - 3, my + 2, vcx + v1x, vcy + v1y, prodHot));
          spawned.push(makeMol(rx.p[2], mx, my - 3, vcx + v2x, vcy + v2y, prodHot));
        }
        const fCol = rx.cat === "ramificación" ? "#ffd43b"
                   : rx.cat === "cierre CO₂" ? "#3fb950"
                   : rx.dH < -150 ? "#ffa657" : "#79c0ff";
        state.flashes.push({ x: mx, y: my, age: 0, col: fCol, r: rx.dH < -150 ? 28 : 18 });
        state.rxnCount++;
        reacted = true;
        appendRxnLog(rx);
        break;
      }

      // Rebote elástico 2D en contacto físico real (rSum = a.r + b.r)
      const rPhys = a.r + b.r;
      if (d2 < rPhys * rPhys) {
        const d = Math.sqrt(d2) || 0.001;
        const nx = dx / d, ny = dy / d;
        const push = (rPhys - d) * 0.5;
        a.x -= nx * push; a.y -= ny * push;
        b.x += nx * push; b.y += ny * push;
        const dvx = b.vx - a.vx, dvy = b.vy - a.vy;
        const vn = dvx * nx + dvy * ny;
        if (vn < 0) {
          const imp = (2 * vn) / (a.M + b.M);
          a.vx += imp * b.M * nx; a.vy += imp * b.M * ny;
          b.vx -= imp * a.M * nx; b.vy -= imp * a.M * ny;
        }
      }
    }
  }

  if (dead.size || spawned.length) {
    state.parts = state.parts.filter(p => !dead.has(p));
    for (const s of spawned) state.parts.push(s);
  }

  if (reacted) {
    const cNow = currentCounts();
    const Told = state.T;
    state.T = computeTemperature(cNow, state.Hchem0, state.sparkKJ, state.lossKJ, state.initO2, state.T0);
    rescaleVelocities(Told, state.T);
  }
}

function stepReactor(dt) {
  // 2 subpasos para captura suave de colisiones rápidas en la cámara 2D
  stepSub(dt * 0.5);
  stepSub(dt * 0.5);

  for (let i = state.flashes.length - 1; i >= 0; i--) {
    state.flashes[i].age += dt;
    if (state.flashes[i].age > 0.55) state.flashes.splice(i, 1);
  }
  if (state.sparkFX) {
    state.sparkFX.age += dt;
    if (state.sparkFX.age > 0.65) state.sparkFX = null;
  }

  // Muestreo temporal uniforme de historia y refresco UI a ~10 Hz (evita thrash DOM)
  histAccum += dt;
  uiAccum += dt;
  if (histAccum >= 0.15) {
    histAccum -= 0.15;
    const cNow = currentCounts();
    recordHistory(cNow);
    drawKinetics();
    drawProfile(cNow);
  }
  if (uiAccum >= 0.10) {
    uiAccum = 0;
    updateSideUI(currentCounts());
  }
  state.reactorDirty = true;
}

// ---------- Dibujo de moléculas en el reactor 2D ----------
function drawAtomBall(ctx, x, y, r, fill) {
  ctx.fillStyle = fill;
  ctx.beginPath(); ctx.arc(x, y, r, 0, TAU); ctx.fill();
  ctx.fillStyle = "rgba(255,255,255,0.32)";
  ctx.beginPath(); ctx.arc(x - r * 0.28, y - r * 0.28, r * 0.34, 0, TAU); ctx.fill();
}

function drawMolGlyph(ctx, p) {
  const s = SPECIES[p.sp];
  ctx.save();
  ctx.translate(p.x, p.y);
  ctx.rotate(p.ang);

  if (s.rad) {
    ctx.strokeStyle = s.col + "aa";
    ctx.lineWidth = 1.4;
    ctx.setLineDash([2.5, 2]);
    ctx.beginPath(); ctx.arc(0, 0, p.r + 3.0, 0, TAU); ctx.stroke();
    ctx.setLineDash([]);
  }

  if (p.sp === "CH4") {
    drawAtomBall(ctx, 0, 0, 5.2, "#6e7681");
    for (let k = 0; k < 4; k++) {
      const a = k * Math.PI * 0.5 + 0.4;
      drawAtomBall(ctx, Math.cos(a) * 6.4, Math.sin(a) * 6.4, 2.7, "#e6edf3");
    }
  } else if (p.sp === "O2") {
    ctx.strokeStyle = "#f8514999"; ctx.lineWidth = 2.6;
    ctx.beginPath(); ctx.moveTo(-4.2, 0); ctx.lineTo(4.2, 0); ctx.stroke();
    drawAtomBall(ctx, -4.4, 0, 4.3, "#f85149");
    drawAtomBall(ctx,  4.4, 0, 4.3, "#f85149");
  } else if (p.sp === "CO2") {
    ctx.strokeStyle = "#3fb95099"; ctx.lineWidth = 2.6;
    ctx.beginPath(); ctx.moveTo(-7.2, 0); ctx.lineTo(7.2, 0); ctx.stroke();
    drawAtomBall(ctx, -7.2, 0, 3.8, "#f85149");
    drawAtomBall(ctx,  7.2, 0, 3.8, "#f85149");
    drawAtomBall(ctx,  0,   0, 4.4, "#3fb950");
  } else if (p.sp === "H2O") {
    drawAtomBall(ctx, 0, -1, 4.6, "#58a6ff");
    drawAtomBall(ctx, -4.5, 3.8, 2.6, "#e6edf3");
    drawAtomBall(ctx,  4.5, 3.8, 2.6, "#e6edf3");
  } else if (p.sp === "CO") {
    ctx.strokeStyle = "#ffa657"; ctx.lineWidth = 3.0;
    ctx.beginPath(); ctx.moveTo(-4, 0); ctx.lineTo(4, 0); ctx.stroke();
    drawAtomBall(ctx, -4.2, 0, 4.2, "#ffa657");
    drawAtomBall(ctx,  4.2, 0, 4.2, "#f85149");
  } else if (p.sp === "CH2O") {
    drawAtomBall(ctx, -1, 0, 4.4, "#c084fc");
    drawAtomBall(ctx,  5.6, 0, 3.8, "#f85149");
    drawAtomBall(ctx, -5.2, -3.5, 2.4, "#e6edf3");
    drawAtomBall(ctx, -5.2,  3.5, 2.4, "#e6edf3");
  } else if (p.sp === "H2") {
    drawAtomBall(ctx, -2.8, 0, 2.8, "#e6edf3");
    drawAtomBall(ctx,  2.8, 0, 2.8, "#e6edf3");
  } else {
    drawAtomBall(ctx, 0, 0, p.r * 0.72, s.col);
    ctx.rotate(-p.ang);
    ctx.fillStyle = "#0d1117";
    ctx.font = "bold 7.5px ui-monospace,monospace";
    ctx.textAlign = "center";
    ctx.fillText(s.label.replace(/·/g, ""), 0, 2.8);
  }
  ctx.restore();
}

function drawReactor() {
  const { W, H } = fitCanvas(cvReactor, ctxR);
  if (!W || !H) return;
  state.W = W; state.H = H;
  ctxR.fillStyle = "#090d13";
  ctxR.fillRect(0, 0, W, H);

  const cNow = currentCounts();
  const nCO = cNow.CO || 0, nCO2 = cNow.CO2 || 0;
  const nRad = (cNow.CH3 || 0) + (cNow.OH || 0) + (cNow.O || 0) + (cNow.H || 0) + (cNow.HO2 || 0) + (cNow.CHO || 0);
  const coRatio = (nCO + nCO2) > 0 ? nCO / (nCO + nCO2) : 0;
  const reg = flameRegime(state.phi, state.T, coRatio, state.rxnCount, cNow.O2 || 0, nRad);
  const g = ctxR.createRadialGradient(W * 0.5, H * 0.5, 16, W * 0.5, H * 0.5, Math.max(W, H) * 0.65);
  g.addColorStop(0, reg.glow);
  g.addColorStop(1, "rgba(9,13,19,0)");
  ctxR.fillStyle = g;
  ctxR.fillRect(0, 0, W, H);

  for (const p of state.parts) drawMolGlyph(ctxR, p);

  for (const f of state.flashes) {
    const u = f.age / 0.55;
    ctxR.strokeStyle = f.col;
    ctxR.globalAlpha = (1 - u) * 0.85;
    ctxR.lineWidth = 2.2 * (1 - u) + 0.5;
    ctxR.beginPath(); ctxR.arc(f.x, f.y, f.r * (0.35 + u), 0, TAU); ctxR.stroke();
    ctxR.globalAlpha = 1;
  }

  if (state.sparkFX) {
    const u = state.sparkFX.age / 0.65;
    const { x, y } = state.sparkFX;
    ctxR.save();
    ctxR.globalAlpha = 1 - u;
    ctxR.strokeStyle = "#fff3bf";
    ctxR.lineWidth = 2.5;
    for (let k = 0; k < 8; k++) {
      const a = k * Math.PI / 4 + u * 0.4;
      const r1 = 6, r2 = 28 + u * 36;
      ctxR.beginPath();
      ctxR.moveTo(x + Math.cos(a) * r1, y + Math.sin(a) * r1);
      ctxR.lineTo(x + Math.cos(a + 0.14) * ((r1 + r2) * 0.5), y + Math.sin(a + 0.14) * ((r1 + r2) * 0.5));
      ctxR.lineTo(x + Math.cos(a) * r2, y + Math.sin(a) * r2);
      ctxR.stroke();
    }
    ctxR.restore();
  }

  // Pastilla de estado de llama usando altura cacheada (cero reflow por frame)
  ctxR.save();
  ctxR.font = "bold 11px system-ui,sans-serif";
  const maxW = W - 24;
  let txt = reg.label;
  while (txt.length > 18 && ctxR.measureText(txt).width + 18 > maxW) txt = txt.slice(0, -4) + "…";
  const tw = Math.min(maxW, ctxR.measureText(txt).width + 18);
  const bx = Math.max(12, W / 2 - tw / 2), by = Math.max(26, state.infoH + 10);
  ctxR.fillStyle = "#0d1117e6";
  ctxR.strokeStyle = reg.col;
  ctxR.lineWidth = 1.4;
  ctxR.beginPath(); roundRectPath(ctxR, bx, by, tw, 22, 11); ctxR.fill(); ctxR.stroke();
  ctxR.fillStyle = reg.col; ctxR.textAlign = "center";
  ctxR.fillText(txt, bx + tw / 2, by + 15);
  ctxR.restore();
}

// ---------- Gráfica 1: Cinética de especies y T(t) ----------
function drawKinetics() {
  const { W, H } = fitCanvas(cvKin, ctxK);
  if (!W || !H) return;
  ctxK.fillStyle = "#0d1117"; ctxK.fillRect(0, 0, W, H);
  const padL = 28, padR = 36, padT = 20, padB = 20;
  const pW = W - padL - padR, pH = H - padT - padB;
  const hist = state.history;
  const maxN = Math.max(state.initO2, state.initCH4 * 2, 20);
  const tMin = hist.length ? hist[0].t : 0;
  const tMax = Math.max(tMin + 4, hist.length ? hist[hist.length - 1].t : 4);
  const tMinK = 250, tMaxK = 2750;
  const X = t => padL + ((t - tMin) / (tMax - tMin)) * pW;
  const Yn = n => padT + (1 - Math.min(1, n / maxN)) * pH;
  const Yt = T => padT + (1 - Math.max(0, Math.min(1, (T - tMinK) / (tMaxK - tMinK)))) * pH;

  ctxK.strokeStyle = "#1a2029"; ctxK.lineWidth = 1;
  for (let i = 1; i <= 3; i++) {
    const y = padT + (i / 4) * pH;
    ctxK.beginPath(); ctxK.moveTo(padL, y); ctxK.lineTo(W - padR, y); ctxK.stroke();
  }

  const series = [
    ["CH4", "#6e7681"], ["O2", "#f85149"], ["CO2", "#3fb950"],
    ["H2O", "#58a6ff"], ["CO", "#ffa657"], ["rad", "#b197fc"]
  ];
  if (hist.length >= 2) {
    for (const [k, col] of series) {
      ctxK.beginPath();
      hist.forEach((pt, i) => {
        const x = X(pt.t), y = Yn(pt[k]);
        i === 0 ? ctxK.moveTo(x, y) : ctxK.lineTo(x, y);
      });
      ctxK.strokeStyle = col;
      ctxK.lineWidth = k === "CO" || k === "CO2" ? 1.9 : 1.35;
      ctxK.stroke();
    }
    ctxK.save(); ctxK.setLineDash([4, 3]);
    ctxK.beginPath();
    hist.forEach((pt, i) => {
      const x = X(pt.t), y = Yt(pt.T);
      i === 0 ? ctxK.moveTo(x, y) : ctxK.lineTo(x, y);
    });
    ctxK.strokeStyle = "#f0c000"; ctxK.lineWidth = 1.6; ctxK.stroke();
    ctxK.restore();
  }

  ctxK.strokeStyle = "#30363d";
  ctxK.beginPath(); ctxK.moveTo(padL, padT); ctxK.lineTo(padL, H - padB); ctxK.lineTo(W - padR, H - padB); ctxK.stroke();
  ctxK.fillStyle = "#8b949e"; ctxK.font = "9px ui-monospace,monospace";
  ctxK.textAlign = "right";
  ctxK.fillText(String(maxN), padL - 3, padT + 4);
  ctxK.fillText("0", padL - 3, H - padB + 3);
  ctxK.textAlign = "left"; ctxK.fillStyle = "#f0c000";
  ctxK.fillText("2750", W - padR + 3, padT + 4);
  ctxK.fillText("250K", W - padR + 3, H - padB + 3);
  ctxK.textAlign = "center"; ctxK.fillStyle = "#8b949e";
  ctxK.fillText(pW < 185 ? "t (sim)" : "Tiempo simulado · especies y T(t)", padL + pW / 2, H - 5);
  ctxK.textAlign = "left"; ctxK.fillStyle = "#e6edf3"; ctxK.font = "bold 10.5px system-ui,sans-serif";
  ctxK.fillText("Cinética y T(t)", padL + 2, 13);
}

// ---------- Gráfica 2: Perfil de entalpía H(ξ) ----------
function drawProfile(cNow) {
  const { W, H } = fitCanvas(cvProf, ctxP);
  if (!W || !H) return;
  ctxP.fillStyle = "#0d1117"; ctxP.fillRect(0, 0, W, H);
  const padL = 38, padR = 12, padT = 22, padB = 22;
  const pW = W - padL - padR, pH = H - padT - padB;

  const hMin = -980, hMax = 380;
  const X = u => padL + u * pW;
  const Y = h => padT + (1 - (h - hMin) / (hMax - hMin)) * pH;

  ctxP.strokeStyle = "#21262d"; ctxP.setLineDash([3, 3]);
  ctxP.beginPath(); ctxP.moveTo(padL, Y(0)); ctxP.lineTo(W - padR, Y(0)); ctxP.stroke();
  ctxP.setLineDash([]);

  const knots = [
    [0.00, 0], [0.14, 0], [0.27, 238], [0.42, -90],
    [0.54, -220], [0.68, -510], [0.80, -802.5], [1.00, -802.5]
  ];
  ctxP.beginPath();
  ctxP.moveTo(X(knots[0][0]), Y(knots[0][1]));
  for (let i = 1; i < knots.length; i++) {
    const [x0, y0] = knots[i - 1], [x1, y1] = knots[i];
    const xm = (x0 + x1) * 0.5;
    ctxP.bezierCurveTo(X(xm), Y(y0), X(xm), Y(y1), X(x1), Y(y1));
  }
  ctxP.strokeStyle = "#58a6ff"; ctxP.lineWidth = 2; ctxP.stroke();

  const xArr = X(0.89);
  ctxP.strokeStyle = "#3fb950"; ctxP.lineWidth = 1.4;
  ctxP.beginPath(); ctxP.moveTo(xArr, Y(0)); ctxP.lineTo(xArr, Y(-802.5)); ctxP.stroke();
  ctxP.fillStyle = "#3fb950";
  ctxP.beginPath(); ctxP.moveTo(xArr, Y(-802.5)); ctxP.lineTo(xArr - 3.5, Y(-802.5) - 6); ctxP.lineTo(xArr + 3.5, Y(-802.5) - 6); ctxP.closePath(); ctxP.fill();
  ctxP.font = "bold 9.5px ui-monospace,monospace"; ctxP.textAlign = "right";
  ctxP.fillText("−802.5", xArr - 5, Y(-420));

  const dHperCH4 = (chemicalEnthalpy(cNow) - state.Hchem0) / Math.max(1, state.initCH4);
  const xi = Math.max(0, Math.min(1, -dHperCH4 / 802.5));
  const curX = X(0.08 + xi * 0.80);
  const curY = Y(Math.max(hMin + 25, Math.min(hMax - 20, dHperCH4)));
  ctxP.fillStyle = "#f0c000";
  ctxP.beginPath(); ctxP.arc(curX, curY, 4.4, 0, TAU); ctxP.fill();
  ctxP.strokeStyle = "#0d1117"; ctxP.lineWidth = 1.2; ctxP.stroke();

  ctxP.fillStyle = "#8b949e"; ctxP.font = "9px ui-monospace,monospace";
  ctxP.textAlign = "center";
  ctxP.fillText("Eₐ=238", X(0.27), Y(238) - 5);
  ctxP.fillText("CH₄+2O₂", X(0.13), Y(0) + 12);
  ctxP.fillText("CO₂+2H₂O", X(0.82), Math.min(H - padB - 4, Y(-802.5) + 12));

  ctxP.strokeStyle = "#30363d";
  ctxP.beginPath(); ctxP.moveTo(padL, padT); ctxP.lineTo(padL, H - padB); ctxP.lineTo(W - padR, H - padB); ctxP.stroke();
  ctxP.textAlign = "right";
  ctxP.fillText("+238", padL - 3, Y(238) + 3);
  ctxP.fillText("0", padL - 3, Y(0) + 3);
  ctxP.fillText("−803", padL - 3, Y(-802.5) + 3);
  ctxP.textAlign = "left"; ctxP.fillStyle = "#e6edf3"; ctxP.font = "bold 10.5px system-ui,sans-serif";
  ctxP.fillText("Entalpía H(ξ) [kJ/mol]", padL + 2, 13);
}

// ---------- Gráfica 3: Visor 3D interactivo (sin solape con barra de botones) ----------
const EL_STYLE = {
  C: { r: 16, col: "#7d8590", hi: "#c9d1d9" },
  O: { r: 15, col: "#f85149", hi: "#ffa198" },
  H: { r: 10, col: "#dce3ea", hi: "#ffffff" }
};

function drawMol3D() {
  const { W, H } = fitCanvas(cvMol, ctxM);
  if (!W || !H) return;
  ctxM.fillStyle = "#0d1117"; ctxM.fillRect(0, 0, W, H);
  const mol = MOL_3D[state.molKey] || MOL_3D.CH4;
  const cx = W * 0.5, cy = H * 0.49;
  const scale = Math.min(W, H) * 0.25;
  const cY = Math.cos(state.rotY), sY = Math.sin(state.rotY);
  const cX = Math.cos(state.rotX), sX = Math.sin(state.rotX);

  const proj = mol.atoms.map((a, idx) => {
    const x1 = a.x * cY + a.z * sY;
    const z1 = -a.x * sY + a.z * cY;
    const y2 = a.y * cX - z1 * sX;
    const z2 = a.y * sX + z1 * cX; // +z2 = más cerca del observador
    const persp = 1 / (1 - z2 * 0.14);
    const st = EL_STYLE[a.el];
    const rad = st.r * persp * Math.min(1.05, Math.min(W, H) / 165);
    return { idx, el: a.el, radical: a.radical, sx: cx + x1 * scale * persp, sy: cy - y2 * scale * persp, z: z2, rad };
  });

  const items = [];
  for (const [i, j, order] of mol.bonds) {
    const a = proj[i], b = proj[j];
    items.push({ type: "bond", z: (a.z + b.z) * 0.5, a, b, order });
  }
  for (const a of proj) items.push({ type: "atom", z: a.z, a });
  items.sort((u, v) => u.z - v.z);

  for (const it of items) {
    if (it.type === "bond") {
      const { a, b, order } = it;
      const dx = b.sx - a.sx, dy = b.sy - a.sy;
      const len = Math.hypot(dx, dy) || 1;
      const ux = dx / len, uy = dy / len;
      const nx = -uy, ny = ux;
      // Recorta ligeramente el enlace para que termine limpio en las esferas
      const x1 = a.sx + ux * a.rad * 0.55, y1 = a.sy + uy * a.rad * 0.55;
      const x2 = b.sx - ux * b.rad * 0.55, y2 = b.sy - uy * b.rad * 0.55;
      const offs = order === 1 ? [0] : order === 2 ? [-3.2, 3.2] : [-4.4, 0, 4.4];
      ctxM.strokeStyle = "#8b949ecc";
      ctxM.lineWidth = order === 1 ? 3.2 : 2.1;
      ctxM.lineCap = "round";
      for (const o of offs) {
        ctxM.beginPath();
        ctxM.moveTo(x1 + nx * o, y1 + ny * o);
        ctxM.lineTo(x2 + nx * o, y2 + ny * o);
        ctxM.stroke();
      }
    } else {
      const { a } = it;
      const st = EL_STYLE[a.el];
      const rad = a.rad;
      const grad = ctxM.createRadialGradient(a.sx - rad * 0.3, a.sy - rad * 0.3, rad * 0.12, a.sx, a.sy, rad);
      grad.addColorStop(0, st.hi);
      grad.addColorStop(0.55, st.col);
      grad.addColorStop(1, "#0d1117");
      ctxM.fillStyle = grad;
      ctxM.beginPath(); ctxM.arc(a.sx, a.sy, rad, 0, TAU); ctxM.fill();
      ctxM.fillStyle = a.el === "H" ? "#0d1117" : "#fff";
      ctxM.font = `bold ${Math.max(8, Math.round(rad * 0.72))}px system-ui,sans-serif`;
      ctxM.textAlign = "center";
      ctxM.fillText(a.el, a.sx, a.sy + rad * 0.25);
      if (a.radical) {
        ctxM.fillStyle = "#f0c000";
        ctxM.beginPath(); ctxM.arc(a.sx - rad * 1.15, a.sy - rad * 0.5, 3.4, 0, TAU); ctxM.fill();
      }
    }
  }

  ctxM.fillStyle = "#e6edf3"; ctxM.font = "bold 10.5px system-ui,sans-serif"; ctxM.textAlign = "left";
  ctxM.fillText(`${mol.formula} — ${mol.name}`, 8, 13);
  ctxM.fillStyle = "#8b949e"; ctxM.font = "9.5px ui-monospace,monospace"; ctxM.textAlign = "center";
  // Dos líneas limpias abajo para nunca recortarse en anchos ~280–305 px
  ctxM.fillText(mol.geom, W * 0.5, H - 16);
  ctxM.fillText(mol.bondInfo, W * 0.5, H - 4);
  state.molDirty = false;
}

// ---------- Panel lateral ----------
function updateSideUI(c) {
  const atoms = countAtoms(c);
  const dHrel = state.Hchem0 - chemicalEnthalpy(c);
  const nCO = c.CO || 0, nCO2 = c.CO2 || 0;
  const nRad = (c.CH3 || 0) + (c.OH || 0) + (c.O || 0) + (c.H || 0) + (c.HO2 || 0) + (c.CHO || 0);
  const Tth = adiabaticFlameTemp(state.initCH4, state.initO2, state.T0);
  const preset = PRESETS.find(pr => pr.id === state.presetId);

  $("rxnDesc").textContent = preset
    ? preset.desc
    : `Mezcla personalizada φ = ${state.phi.toFixed(2)} (${state.initCH4} CH₄ + ${state.initO2} O₂ en aire).`;
  $("phiVal").textContent = `φ = ${state.phi.toFixed(2)}`;
  $("t0Val").textContent = `${state.T0} K`;

  $("stCH4").textContent = c.CH4 || 0;
  $("stCO2").textContent = nCO2;
  $("stCO").textContent = nCO;
  $("stRad").textContent = nRad;

  $("kvT").textContent = `${Math.round(state.T)} K (${Math.round(state.T - 273.15)} °C)`;
  $("kvTad").textContent = `~${Math.round(Tth)} K`;
  const dhStr = Math.abs(dHrel) < 0.5 ? "0 kJ" : `${dHrel > 0 ? "−" : "+"}${Math.abs(dHrel).toFixed(0)} kJ`;
  $("kvDH").textContent = `${dhStr} (en ${state.initCH4} mol eq.)`;
  $("kvAtoms").textContent = `C:${atoms.C}/${state.initAtoms.C} · H:${atoms.H}/${state.initAtoms.H} · O:${atoms.O}/${state.initAtoms.O}`;
  $("kvRxns").textContent = String(state.rxnCount);

  const infoEl = $("reactorInfo");
  infoEl.textContent =
    `T = ${Math.round(state.T)} K · φ = ${state.phi.toFixed(2)} · CH₄:${c.CH4 || 0} O₂:${c.O2 || 0} → CO₂:${nCO2} H₂O:${c.H2O || 0} CO:${nCO}`;
  state.infoH = infoEl.offsetHeight || 18;
}

function setPaused(p) {
  state.paused = !!p;
  const b = $("btnPause");
  b.textContent = state.paused ? "▶ Reanudar" : "⏸ Pausa";
  b.setAttribute("aria-pressed", String(state.paused));
}

function selectPreset(id) {
  const pr = PRESETS.find(p => p.id === id) || PRESETS[0];
  state.presetId = pr.id;
  state.phi = pr.phi;
  state.T0 = pr.T0;
  $("phiSlider").value = String(pr.phi);
  $("t0Slider").value = String(pr.T0);
  document.querySelectorAll("#tabs .tab").forEach(el => {
    const on = el.dataset.id === pr.id;
    el.classList.toggle("on", on);
    el.setAttribute("aria-pressed", String(on));
  });
  setPaused(false);
  resetReactor();
}

function buildUI() {
  const tabs = $("tabs");
  PRESETS.forEach(pr => {
    const b = document.createElement("button");
    b.className = "tab";
    b.dataset.id = pr.id;
    const strong = document.createElement("b");
    strong.textContent = pr.label;
    const sm = document.createElement("small");
    sm.textContent = pr.sub;
    b.append(strong, sm);
    b.onclick = () => selectPreset(pr.id);
    tabs.appendChild(b);
  });

  const molBar = $("molBar");
  Object.keys(MOL_3D).forEach(k => {
    const b = document.createElement("button");
    b.className = "mbtn" + (k === state.molKey ? " on" : "");
    b.dataset.m = k;
    b.setAttribute("aria-pressed", String(k === state.molKey));
    b.textContent = MOL_3D[k].formula;
    b.onclick = () => {
      state.molKey = k;
      document.querySelectorAll(".mbtn").forEach(x => {
        const on = x.dataset.m === k;
        x.classList.toggle("on", on);
        x.setAttribute("aria-pressed", String(on));
      });
      state.molDirty = true;
      drawMol3D();
    };
    molBar.appendChild(b);
  });

  // Actualiza lectura en oninput; reinicia reactor tras breve debounce o en onchange
  let sliderTimer = null;
  const scheduleReset = () => {
    clearTimeout(sliderTimer);
    sliderTimer = setTimeout(resetReactor, 110);
  };
  $("phiSlider").oninput = e => {
    state.presetId = "custom";
    state.phi = parseFloat(e.target.value);
    $("phiVal").textContent = `φ = ${state.phi.toFixed(2)}`;
    document.querySelectorAll("#tabs .tab").forEach(x => { x.classList.remove("on"); x.setAttribute("aria-pressed", "false"); });
    scheduleReset();
  };
  $("phiSlider").onchange = () => { clearTimeout(sliderTimer); resetReactor(); };

  $("t0Slider").oninput = e => {
    state.presetId = "custom";
    state.T0 = parseInt(e.target.value, 10);
    $("t0Val").textContent = `${state.T0} K`;
    document.querySelectorAll("#tabs .tab").forEach(x => { x.classList.remove("on"); x.setAttribute("aria-pressed", "false"); });
    scheduleReset();
  };
  $("t0Slider").onchange = () => { clearTimeout(sliderTimer); resetReactor(); };

  $("chkAdiab").onchange = e => { state.adiabatic = e.target.checked; };
  $("btnSpark").onclick = fireSpark;
  $("btnPause").onclick = () => setPaused(!state.paused);
  $("btnReset").onclick = () => { setPaused(false); resetReactor(); };

  document.querySelectorAll(".spd").forEach(b => {
    const init = parseFloat(b.dataset.s) === state.speed;
    b.classList.toggle("on", init);
    b.setAttribute("aria-pressed", String(init));
    b.onclick = () => {
      document.querySelectorAll(".spd").forEach(x => { x.classList.remove("on"); x.setAttribute("aria-pressed", "false"); });
      b.classList.add("on"); b.setAttribute("aria-pressed", "true");
      state.speed = parseFloat(b.dataset.s);
    };
  });

  // Arrastrar y teclado en el visor 3D
  let drag3d = null;
  cvMol.addEventListener("pointerdown", e => {
    drag3d = { x: e.clientX, y: e.clientY };
    cvMol.setPointerCapture(e.pointerId);
  });
  cvMol.addEventListener("pointermove", e => {
    if (!drag3d) return;
    state.rotY += (e.clientX - drag3d.x) * 0.014;
    state.rotX = Math.max(-1.4, Math.min(1.4, state.rotX + (e.clientY - drag3d.y) * 0.014));
    drag3d = { x: e.clientX, y: e.clientY };
    state.molDirty = true;
    drawMol3D();
  });
  const endDrag = () => { drag3d = null; };
  cvMol.addEventListener("pointerup", endDrag);
  cvMol.addEventListener("pointercancel", endDrag);
  cvMol.addEventListener("keydown", e => {
    if (e.key === "ArrowLeft")  { e.preventDefault(); state.rotY -= 0.18; state.molDirty = true; drawMol3D(); }
    if (e.key === "ArrowRight") { e.preventDefault(); state.rotY += 0.18; state.molDirty = true; drawMol3D(); }
    if (e.key === "ArrowUp")    { e.preventDefault(); state.rotX = Math.max(-1.4, state.rotX - 0.15); state.molDirty = true; drawMol3D(); }
    if (e.key === "ArrowDown")  { e.preventDefault(); state.rotX = Math.min(1.4, state.rotX + 0.15); state.molDirty = true; drawMol3D(); }
  });

  window.addEventListener("resize", () => {
    state.infoH = $("reactorInfo").offsetHeight || 18;
    state.reactorDirty = true;
    state.molDirty = true;
    drawReactor();
    drawKinetics();
    drawProfile(currentCounts());
    drawMol3D();
  });
}

let last = performance.now();
function loop(now) {
  const dt = Math.min(0.05, (now - last) / 1000) * state.speed;
  last = now;
  if (!document.hidden && !state.paused) {
    stepReactor(dt);
    if (!reducedMotion) {
      state.rotY += dt * 0.35;
      state.molDirty = true;
    }
  }
  if (!document.hidden && (state.reactorDirty || !state.paused)) {
    drawReactor();
    state.reactorDirty = false;
  }
  if (!document.hidden && state.molDirty) {
    drawMol3D();
  }
  requestAnimationFrame(loop);
}

buildUI();
selectPreset("stoich");
requestAnimationFrame(loop);
})();

// ==================== EXPORTS PARA PRUEBAS EN NODE ====================
if (typeof module !== "undefined") {
  module.exports = {
    R_GAS, N2_PER_O2, CP_N2, DH_VAP_H2O,
    SPECIES, BOND_ENERGIES, REACTIONS, UNIMOL_CHO, MOL_3D, PRESETS,
    bondBalance, hessReactionEnthalpy, arrheniusProb,
    equivalenceRatio, mixtureCounts, countAtoms,
    chemicalEnthalpy, totalHeatCapacity, computeTemperature,
    adiabaticFlameTemp, flameRegime
  };
}
})();
