#!/usr/bin/env bash
# Baja a Bibliografia/pdf/ los PDFs de arXiv de la bibliografía (Bibliografia/README.md) que falten.
# Los PDFs no se versionan (el repositorio es público). Uso: Bibliografia/descargar.sh
set -euo pipefail
DIR=$(cd "$(dirname "$0")" && pwd)/pdf
mkdir -p "$DIR"
while read -r clave arxiv; do
  [ -f "$DIR/$clave.pdf" ] && continue
  echo "Bajando $clave (arXiv:$arxiv)"
  curl -sSL -A "Mozilla/5.0" -o "$DIR/$clave.pdf" "https://arxiv.org/pdf/$arxiv" || echo "  falló $clave"
  sleep 3  # pausa entre descargas, como pide arXiv
done <<'EOF'
Kofman1997 hep-ph/9704452
Greene1997 hep-ph/9705347
Amin2015 1410.3808
Figueroa2017a 1609.05197
Mishra2024 2403.10606
Barman2024 2404.16090
Elia2025 2510.08685
Podolsky2006 hep-ph/0507096
Lozanov2017 1608.01213
Dufaux2007 0707.0875
GarciaBellido2007 astro-ph/0701014
Figueroa2017b 1707.04533
Gross2024 2411.04189
Das2025 2508.07442
Figueroa2021 2006.15122
Figueroa2023 2102.01031
Figueroa2024 2312.15056
BaezaBallesteros2025 2512.15627
BaezaBallesteros2026 2607.24978
Ellis2026 2510.18656
Saini2024 2409.05615
Ketov2025 2501.06451
Baumann2009 0907.5424
Caprini2018 1801.04268
Christensen2019 1811.08797
Renzini2022 2202.00178
Watanabe2006 astro-ph/0604176
Ringwald2022 2203.00621
Moore2015 1408.0740
Aggarwal2021 2011.12414
Domcke2023 2306.04496
Berlin2022 2112.11465
Herman2023 2203.15668
Gatti2024 2403.18610
Mohanty2026 2503.06858
Amaral2026 2603.24645
Abac2026 2503.12263
EOF
