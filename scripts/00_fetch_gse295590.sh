#!/usr/bin/env bash
# Phase 2 — download the GSE295590 processed matrices from GEO.
#
#   bash scripts/00_fetch_gse295590.sh            # -> data/raw/GSE295590/
#   OUT_DIR=/scratch/ilc2 bash scripts/00_fetch_gse295590.sh
#
# Downloads GSE295590_RAW.tar (per-sample matrix.mtx.gz / barcodes.tsv.gz /
# features.tsv.gz), records its SHA-256 and the download date in
# download_log.tsv (pipeline Phase 2: "record accession numbers and download
# information"), then extracts it. Raw FASTQ is NOT downloaded: processed
# matrices are sufficient until a step requires raw-level reprocessing.
set -euo pipefail

SERIES="GSE295590"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/${SERIES}/suppl/${SERIES}_RAW.tar"
OUT_DIR="${OUT_DIR:-data/raw/${SERIES}}"
TAR="${OUT_DIR}/${SERIES}_RAW.tar"
LOG="${OUT_DIR}/download_log.tsv"

mkdir -p "${OUT_DIR}"

if [[ -s "${TAR}" ]]; then
  echo "[skip] ${TAR} already present."
else
  echo "[get ] ${URL}"
  curl --fail --location --retry 3 --continue-at - --output "${TAR}.part" "${URL}"
  mv "${TAR}.part" "${TAR}"
fi

# Refuse to continue on a corrupt archive.
tar -tf "${TAR}" > /dev/null || { echo "[fail] ${TAR} is not a valid tar archive." >&2; exit 1; }

if command -v sha256sum > /dev/null; then
  SHA=$(sha256sum "${TAR}" | cut -d' ' -f1)
else
  SHA=$(shasum -a 256 "${TAR}" | cut -d' ' -f1)   # macOS
fi

[[ -f "${LOG}" ]] || printf 'series\tfile\turl\tsha256\tdownloaded_utc\n' > "${LOG}"
printf '%s\t%s\t%s\t%s\t%s\n' "${SERIES}" "$(basename "${TAR}")" "${URL}" "${SHA}" \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "${LOG}"

tar -xf "${TAR}" -C "${OUT_DIR}"

echo "[done] Files in ${OUT_DIR}:"
ls -1 "${OUT_DIR}" | grep -E '^GSM' || echo "  (no GSM-prefixed files found — inspect the archive manually)"
echo "[note] Expect three GSM8953514/15/16 sample sets. Cell counts are checked in Phase 3, not here."
