#!/usr/bin/env python3
"""Resuelve el conflicto de provenance/ al mergear main en la rama de un modelo.

main y cada rama (monomial, t-model, e-model) agregan entradas al final de provenance/numbers.json
y provenance/claims.yaml, así que al hacer `git merge main` en una rama git suele marcar conflicto
en esos dos archivos. La resolución correcta es la unión: lo de la rama más lo nuevo de main.

Uso, desde la raíz del repositorio y con el merge en conflicto:
    git checkout t-model
    git merge main                          # CONFLICT en provenance/...
    python3 code/resolver_procedencia.py
    git add provenance && git commit --no-edit

numbers.json: conserva las claves de la rama y agrega las que sólo están en main; si una clave está
en los dos lados con valores distintos, se detiene (eso no es un conflicto de formato y hay que
mirarlo a mano). claims.yaml: conserva el texto de la rama y agrega al final las afirmaciones de
main cuyo id no está en la rama. Al final comprueba que no haya ids repetidos y que todo número
citado por una afirmación exista.
"""

import json
import subprocess
import sys

import yaml


def show(etapa, ruta):
    """Contenido de un archivo en una etapa del índice durante un merge (:2 rama, :3 main)."""
    r = subprocess.run(["git", "show", f"{etapa}:{ruta}"], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"No hay conflicto de merge en {ruta} (¿corriste git merge main?): {r.stderr.strip()}")
    return r.stdout


def main():
    rama = json.loads(show(":2", "provenance/numbers.json"))
    otra = json.loads(show(":3", "provenance/numbers.json"))
    for k, v in otra.items():
        if k not in rama:
            rama[k] = v
        elif rama[k] != v:
            sys.exit(f"numbers.json[{k}] tiene valores distintos en la rama y en main: resolvelo a mano.")
    with open("provenance/numbers.json", "w", encoding="utf-8") as fh:
        json.dump(rama, fh, indent=2, ensure_ascii=False)

    txt_rama, txt_otra = show(":2", "provenance/claims.yaml"), show(":3", "provenance/claims.yaml")
    ids_rama = {c["id"] for c in yaml.safe_load(txt_rama)["claims"]}
    bloques = ("\n" + txt_otra).split("\n  - id: ")[1:]
    nuevos = ["  - id: " + b.rstrip("\n") + "\n" for b in bloques if b.split("\n")[0].strip() not in ids_rama]
    txt = txt_rama.rstrip("\n") + "\n" + "".join(nuevos)
    with open("provenance/claims.yaml", "w", encoding="utf-8") as fh:
        fh.write(txt)

    ids = [c["id"] for c in yaml.safe_load(txt)["claims"]]
    faltan = [k for c in yaml.safe_load(txt)["claims"] for k in c.get("numbers", []) if k not in rama]
    if len(ids) != len(set(ids)) or faltan:
        sys.exit(f"Resultado inconsistente: ids repetidos o números faltantes {faltan}")
    print(f"provenance/ resuelto: {len(ids)} afirmaciones ({len(nuevos)} nuevas de main), "
          f"{len(rama)} números. Ahora: git add provenance && git commit --no-edit")


if __name__ == "__main__":
    main()
