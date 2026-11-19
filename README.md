# molecular-simulation

Simulaciones interactivas de química orgánica y reactividad molecular (HTML + CSS + JavaScript puro, sin dependencias ni *build*). Repositorio hermano de [atomic-simulation](https://github.com/alrobles/atomic-simulation).

## Vistas

| Página | Qué muestra |
|---|---|
| [`index.html`](index.html) | **Portada**: modelo 3D interactivo del metano tetraédrico (sp³, 109.5°), catálogo de simulaciones orgánicas, panel de comprobaciones y puente directo hacia `atomic-simulation`. |
| [`combustion.html`](combustion.html) | **Combustión del metano (`CH₄ + 2 O₂ → CO₂ + 2 H₂O`)**: reactor 2D con mecanismo elemental reducido e irreversible en cadena ramificada (29 reacciones bimoleculares + 1 disociación térmica de `·CHO` = 30 pasos con `Eₐ ≥ ΔH°`), chispa eléctrica local de ignición, control de relación de equivalencia `φ` (llama azul completa vs combustión rica incompleta con `CO` tóxico y `H₂`), Primera Ley por Ley de Hess (`ΔH° = −802.5 kJ/mol` en gas / `−890.5 kJ/mol` líquido; `T_ad ≈ 2235 K` en aire con `C_p` medio de llama sin disociación de equilibrio), gráficas de cinética y entalpía `H(ξ)` y visor 3D interactivo de geometrías moleculares (`CH₄`, `O₂`, `CH₂O`, `CO`, `CO₂`, `H₂O`, `·OH`). |

## Estructura

```
index.html          portada orgánica con CH₄ 3D interactivo
combustion.html     simulador de combustión del metano
css/tokens.css      tokens de diseño compartidos (familia atomic-simulation)
css/base.css        primitivas de layout, navegación y accesibilidad
js/combustion.js    termoquímica, cinética de radicales, reactor 2D y visor 3D
                    (funciones puras exportadas para pruebas en Node)
tests/check.js      pruebas automatizadas de física y química (npm test)
```

## Pruebas y ejecución local

```bash
npm test            # ejecuta node tests/check.js
npm run serve       # servidor local en http://localhost:8080/
```

El flujo de GitHub Actions (`.github/workflows/test.yml`) verifica la sintaxis de todos los módulos JS y ejecuta `npm test` en cada *push*.

## Publicación (GitHub Pages)

Publicable directamente desde la rama `main`, raíz del repositorio (*Settings → Pages → Deploy from branch*) en `https://alrobles.github.io/molecular-simulation/`.

## Licencia

MIT — ver [`LICENSE`](LICENSE). © 2026 Ángel Luis Robles Fernández.
