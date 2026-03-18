# Megaliths of a Tilted World - V2

**An observational data note. Not peer-reviewed.**

88 ancient stone construction sites across six continents, classified using standard archaeological masonry typology, tested against 35+ physical property filters. The finding: a great circle tilted approximately 30-34 degrees from the current equator fits the distribution of ancient stone sites better than any other orientation. The pattern is preserved by engineering-quality filters and disrupted by mass, geography, and altitude filters. It tightens as engineering quality increases.

## Quick Start

```bash
python analysis.py              # Run all (2° grid, matches report exactly)
python analysis.py --fast       # Run all (4° grid, ~4x faster, slightly different numbers)
python analysis.py --section 4  # Run only Section 4 (filters)
python analysis.py --quick      # Skip Monte Carlo
```

Requires: Python 3.8+, NumPy, matplotlib (for charts only)

### Grid resolution

The report uses a 2° pole grid (16,380 positions). This is the default. The `--fast` flag uses a 4° grid which runs faster but may produce slightly different optimal tilts (within 2-4°) and site counts (within 1-2) due to coarser sampling. For exact reproduction of report numbers, use the default.

## Files

- `sites_88.json` - Complete dataset: 88 sites with coordinates, masonry classification (0-4), precision score (0-4), max block mass, Mohs hardness, region
- `analysis.py` - All calculations from the report (tilt scans, filter tests, Monte Carlo, cluster analysis)
- `charts/` - Full-resolution versions of all charts in the report (200 DPI PNG)

## Dataset Format

Each site in `sites_88.json`:
```json
{
  "name": "Baalbek",
  "lat": 34.007,
  "lon": 36.204,
  "region": "Levant",
  "masonry": 4,
  "precision": 3,
  "max_tons": 1000,
  "mohs": 3,
  "rock_cut": false,
  "notes": "Trilithon: three 800-ton blocks..."
}
```

**Masonry** (standard archaeological classification):
- 0: Not stone
- 1: Rubble
- 2: Cyclopean
- 3: Polygonal
- 4: Ashlar

**Precision** (author-assigned, definitions in report Appendix C):
- 0: None
- 1: Minimal
- 2: Fitted
- 3: Tight
- 4: Sub-millimeter

## Key Results

| Metric | All 88 | Score >= 7 (21 sites) |
|--------|--------|----------------------|
| Optimal tilt | 32 degrees | 28 degrees |
| Closer to tilted | 73% | 95% |
| Median distance | 939 km | 328 km |
| Monte Carlo Z | 7.68 | 7.59 |
| Cluster centroids | 71% | 100% (9/9) |

## Reproducibility

Every number in the report can be reproduced by running `analysis.py`. To test a different dataset: modify `sites_88.json` and rerun. To add a site: append to the JSON array with the same fields.

## License

Data and code: CC BY 4.0. Cite as: Fralse, T. (2026). Megaliths of a Tilted World: How Engineering Quality Reveals a Great Circle Pattern. https://github.com/torbjorn-fralse/megaliths
