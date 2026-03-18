#!/usr/bin/env python3
"""
Megaliths of a Tilted World - Analysis Code
================================================
Reproduces all findings in the report.

Usage:
    python analysis.py              # Run everything (2 deg grid, matches report exactly)
    python analysis.py --fast       # Use 4 deg grid (faster, slightly different numbers)
    python analysis.py --section 1  # Run specific section
    python analysis.py --quick      # Skip Monte Carlo

Grid resolution: The report uses a 2 deg grid (16,380 pole positions). This is
the default. The --fast flag uses a 4 deg grid which runs ~4x faster but may
produce slightly different optimal tilts (typically within 2-4 deg) and site
counts (within 1-2) due to coarser sampling.

Requires: numpy, matplotlib (for charts only)
"""

import json, math, os, sys
import numpy as np

EARTH_R = 6371.0
TDEG500 = math.degrees(500 / EARTH_R)

# Grid resolution: 2 deg matches report, 4 deg is faster
GRID = 2  # default, overridden by --fast

# ============================================================
# Core functions
# ============================================================

def ll2xyz(lat, lon):
    """Convert lat/lon to unit xyz on sphere."""
    la, lo = math.radians(lat), math.radians(lon)
    return np.array([math.cos(la)*math.cos(lo), math.cos(la)*math.sin(lo), math.sin(la)])

def dist_to_gc(lat, lon, pole):
    """Distance (km) from point to great circle defined by pole."""
    xyz = ll2xyz(lat, lon)
    dot = np.dot(xyz, pole)
    return abs(90.0 - math.degrees(math.acos(np.clip(abs(dot), 0, 1)))) * (math.pi/180) * EARTH_R

def point_dist(lat1, lon1, lat2, lon2):
    """Great circle distance between two points (km)."""
    xyz1 = ll2xyz(lat1, lon1)
    xyz2 = ll2xyz(lat2, lon2)
    return math.acos(np.clip(np.dot(xyz1, xyz2), -1, 1)) * EARTH_R

def optimize_median(site_list):
    """Find pole that minimizes median distance to sites."""
    if len(site_list) < 3:
        return None, None, None
    best_med = 99999
    best_pole = None
    for lat in range(-90, 91, GRID):
        for lon in range(-180, 180, GRID):
            p = ll2xyz(lat, lon)
            dts = np.array([dist_to_gc(s["lat"], s["lon"], p) for s in site_list])
            med = float(np.median(dts))
            if med < best_med:
                best_med = med
                best_pole = (lat, lon)
    tilt = 90 - abs(best_pole[0])
    return tilt, best_med, best_pole

def optimize_500(site_list):
    """Find pole that maximizes sites within 500 km."""
    if len(site_list) < 3:
        return None, None, None
    xyzs = np.array([ll2xyz(s["lat"], s["lon"]) for s in site_list])
    best_count = 0
    best_pole = None
    for lat in range(-90, 91, GRID):
        for lon in range(-180, 180, GRID):
            p = ll2xyz(lat, lon)
            dots = np.abs(xyzs @ p)
            angs = np.arccos(np.clip(dots, 0, 1)) * 180 / np.pi
            c = int(np.sum(np.abs(90 - angs) <= TDEG500))
            if c > best_count:
                best_count = c
                best_pole = (lat, lon)
    tilt = 90 - abs(best_pole[0])
    return tilt, best_count, best_pole

def find_clusters(site_list, radius_km=500):
    """Agglomerative clustering by point-to-point distance."""
    n = len(site_list)
    visited = [False] * n
    clusters = []
    for i in range(n):
        if visited[i]:
            continue
        cluster = [i]
        visited[i] = True
        queue = [i]
        while queue:
            curr = queue.pop(0)
            for j in range(n):
                if visited[j]:
                    continue
                d = point_dist(site_list[curr]["lat"], site_list[curr]["lon"],
                              site_list[j]["lat"], site_list[j]["lon"])
                if d <= radius_km:
                    visited[j] = True
                    cluster.append(j)
                    queue.append(j)
        clusters.append(cluster)
    return clusters

# ============================================================
# Load data
# ============================================================

def load_sites(path="sites_88.json"):
    sites = json.load(open(path))
    pole = ll2xyz(-56, 30)
    pole_n = ll2xyz(90, 0)
    for s in sites:
        s["dt"] = dist_to_gc(s["lat"], s["lon"], pole)
        s["dc"] = dist_to_gc(s["lat"], s["lon"], pole_n)
        s["combined"] = s["masonry"] + s["precision"]
    return sites

IN_SITU = {"Yangshan Quarry", "Unfinished Obelisk (Aswan)", "Gommateshwara"}

# ============================================================
# Section 1: Broad dataset tilt scan
# ============================================================

def section1(sites):
    print("\n" + "="*60)
    print("SECTION 1: THE BROAD DATASET")
    print("="*60)
    
    xyzs = np.array([ll2xyz(s["lat"], s["lon"]) for s in sites])
    pole_n = ll2xyz(90, 0)
    
    print(f"\nTilt scan (all 88 sites, 500 km, {GRID} deg grid):")
    print(f"{'Tilt':>5s} {'<500km':>7s} {'<1000km':>8s} {'Closer%':>8s} {'Median':>8s}")
    
    for tilt in range(0, 91, 5):
        plat = 90 - tilt
        best_c = 0
        best_pole = None
        for plon in range(-180, 180, GRID):
            for sign in [1, -1]:
                p = ll2xyz(sign * plat, plon)
                dots = np.abs(xyzs @ p)
                angs = np.arccos(np.clip(dots, 0, 1)) * 180 / np.pi
                c = int(np.sum(np.abs(90 - angs) <= TDEG500))
                if c > best_c:
                    best_c = c
                    best_pole = ll2xyz(sign * plat, plon)
        
        dts = [dist_to_gc(s["lat"], s["lon"], best_pole) for s in sites]
        dcs = [s["dc"] for s in sites]
        w1000 = sum(1 for d in dts if d <= 1000)
        closer = sum(1 for dt, dc in zip(dts, dcs) if dt < dc)
        med = np.median(dts)
        
        print(f"  {tilt:3d}   {best_c:5d}   {w1000:6d}   {closer/len(sites)*100:6.1f}%  {med:7.0f}")

# ============================================================
# Section 2: Identifying the diagonal
# ============================================================

def section2(sites):
    print("\n" + "="*60)
    print("SECTION 2: IDENTIFYING THE DIAGONAL")
    print("="*60)
    
    print("\nOptimized poles by subset:")
    subsets = {
        "All 88": sites,
        "M>=1": [s for s in sites if s["masonry"] >= 1],
        "M>=2": [s for s in sites if s["masonry"] >= 2],
        "M>=3": [s for s in sites if s["masonry"] >= 3],
    }
    
    for name, sl in subsets.items():
        tilt_m, med_m, pole_m = optimize_median(sl)
        tilt_5, count_5, pole_5 = optimize_500(sl)
        print(f"  {name:10s}: median tilt={tilt_m}, 500km tilt={tilt_5}, median={med_m:.0f}km, <500={count_5}")

# ============================================================
# Section 3: Masonry classification
# ============================================================

def section3(sites):
    print("\n" + "="*60)
    print("SECTION 3: MASONRY CLASSIFICATION")
    print("="*60)
    
    for m in range(5):
        sl = [s for s in sites if s["masonry"] >= m]
        cl = sum(1 for s in sl if s["dt"] < s["dc"])
        tilt, med, _ = optimize_median(sl)
        names = ["All (M0+)", "Rubble+ (M1+)", "Cyclopean+ (M2+)", "Polygonal+ (M3+)", "Ashlar (M4)"]
        print(f"  {names[m]:20s}: n={len(sl):2d}, tilt={tilt}, med={med:.0f}km, closer={cl}/{len(sl)} ({cl/len(sl)*100:.0f}%)")
    
    # M4 precision split
    print("\n  M4 precision split:")
    m4_lo = [s for s in sites if s["masonry"] == 4 and s["precision"] < 3]
    m4_hi = [s for s in sites if s["masonry"] == 4 and s["precision"] >= 3]
    cl_lo = sum(1 for s in m4_lo if s["dt"] < s["dc"])
    cl_hi = sum(1 for s in m4_hi if s["dt"] < s["dc"])
    print(f"    M4+P<3:  {cl_lo}/{len(m4_lo)} closer ({cl_lo/len(m4_lo)*100:.0f}%)")
    print(f"    M4+P>=3: {cl_hi}/{len(m4_hi)} closer ({cl_hi/len(m4_hi)*100:.0f}%)")

# ============================================================
# Section 4: Testing every filter
# ============================================================

def section4(sites):
    print("\n" + "="*60)
    print("SECTION 4: TESTING EVERY FILTER")
    print("="*60)
    
    def report(label, sl):
        if len(sl) < 3:
            print(f"  {label:45s} {len(sl):3d}  -- too few --")
            return
        tilt, med, _ = optimize_median(sl)
        cl = sum(1 for s in sl if s["dt"] < s["dc"])
        in_range = "YES" if 28 <= tilt <= 38 else "no"
        print(f"  {label:45s} {len(sl):3d}  tilt={tilt:2d}  med={med:5.0f}  cl={cl/len(sl)*100:5.1f}%  {in_range}")
    
    print(f"\n{'Filter':45s} {'N':>3s}  {'Tilt':>6s}  {'Med':>5s}  {'Cl%':>6s}  {'28-38?':>6s}")
    print("-" * 80)
    
    for m in range(5):
        report(f"M>={m}", [s for s in sites if s["masonry"] >= m])
    for mt in [10, 50, 100, 200]:
        report(f"Mass >= {mt}t", [s for s in sites if s["max_tons"] >= mt and s["name"] not in IN_SITU])
    for mh in [4, 5, 6, 7]:
        report(f"Mohs >= {mh}", [s for s in sites if s["mohs"] >= mh])
    for pr in [2, 3, 4]:
        report(f"Precision >= {pr}", [s for s in sites if s["precision"] >= pr])
    for sc in [5, 6, 7]:
        report(f"Score >= {sc}", [s for s in sites if s["combined"] >= sc])
    for rname, rfilt in [
        ("Europe", lambda s: s["region"] == "Europe"),
        ("Egypt", lambda s: s["region"] == "Egypt"),
        ("Peru+Bolivia", lambda s: s["region"] in {"Peru", "Bolivia"}),
        ("Mesoamerica", lambda s: s["region"] == "Mesoamerica"),
        ("India", lambda s: s["region"] == "India"),
    ]:
        report(rname, [s for s in sites if rfilt(s)])
    report("Mass>=50t + Polygonal", [s for s in sites if s["max_tons"]>=50 and s["masonry"]>=3 and s["name"] not in IN_SITU])
    report("Mass>=50t + Prec>=3", [s for s in sites if s["max_tons"]>=50 and s["precision"]>=3 and s["name"] not in IN_SITU])

# ============================================================
# Section 5: Combined engineering score
# ============================================================

def section5(sites):
    print("\n" + "="*60)
    print("SECTION 5: COMBINED ENGINEERING SCORE")
    print("="*60)
    
    print(f"\n{'Score':>6s} {'N':>4s} {'Tilt':>5s} {'Median':>7s} {'Closer':>10s} {'%':>5s}")
    for sc in range(8):
        sl = [s for s in sites if s["combined"] >= sc]
        if len(sl) < 2: continue
        tilt, med, _ = optimize_median(sl)
        cl = sum(1 for s in sl if s["dt"] < s["dc"])
        print(f"  >={sc:1d}  {len(sl):4d}  {tilt:4d}  {med:6.0f}  {cl:3d}/{len(sl):2d}  {cl/len(sl)*100:5.1f}%")
    
    s7 = [s for s in sites if s["combined"] >= 7]
    print(f"\nScore >= 7 sites ({len(s7)}):")
    for s in sorted(s7, key=lambda x: x["dt"]):
        cl = "T" if s["dt"] < s["dc"] else "C"
        print(f"  {cl} {s['dt']:7.0f}km  M{s['masonry']}+P{s['precision']}={s['combined']}  {s['name']:30s}  {s['region']}")
    
    print("\nPrecision score randomization (10,000 trials)...")
    np.random.seed(42)
    actual_pct = sum(1 for s in s7 if s["dt"] < s["dc"]) / len(s7) * 100
    masonry_vals = [s["masonry"] for s in sites]
    precision_vals = [s["precision"] for s in sites]
    match_count = 0
    for trial in range(10000):
        shuffled_p = np.random.permutation(precision_vals)
        combined = [m + p for m, p in zip(masonry_vals, shuffled_p)]
        high = [i for i, c in enumerate(combined) if c >= 7]
        if len(high) > 0:
            cl = sum(1 for i in high if sites[i]["dt"] < sites[i]["dc"])
            if cl / len(high) * 100 >= actual_pct:
                match_count += 1
    print(f"  Actual: {actual_pct:.1f}%")
    print(f"  Random achieving same: {match_count}/10000 ({match_count/100:.2f}%)")

# ============================================================
# Section 6: Cluster analysis
# ============================================================

def section6(sites):
    print("\n" + "="*60)
    print("SECTION 6: CLUSTER ANALYSIS")
    print("="*60)
    
    s7 = [s for s in sites if s["combined"] >= 7]
    
    print(f"\nClustering at multiple radii (score >= 7):")
    for radius in [100, 200, 300, 500, 750, 1000]:
        cls = find_clusters(s7, radius)
        centroids = []
        for cl in cls:
            lat_m = np.mean([s7[j]["lat"] for j in cl])
            lon_m = np.mean([s7[j]["lon"] for j in cl])
            dt = dist_to_gc(lat_m, lon_m, ll2xyz(-56, 30))
            dc = dist_to_gc(lat_m, lon_m, ll2xyz(90, 0))
            centroids.append({"lat": lat_m, "lon": lon_m, "dt": dt, "dc": dc})
        ct = sum(1 for c in centroids if c["dt"] < c["dc"])
        tilt, med, _ = optimize_median(centroids)
        print(f"  r={radius:5d}km: {len(cls):2d} clusters, {ct}/{len(cls)} closer ({ct/len(cls)*100:.0f}%), tilt={tilt}, med={med:.0f}km")
    
    print(f"\nLeave-one-region-out (score >= 7):")
    from collections import Counter
    regions = Counter(s["region"] for s in s7)
    for region in sorted(regions.keys()):
        remaining = [s for s in s7 if s["region"] != region]
        cl = sum(1 for s in remaining if s["dt"] < s["dc"])
        tilt, med, _ = optimize_median(remaining)
        print(f"  Remove {region:15s} ({regions[region]}): {cl}/{len(remaining)} closer ({cl/len(remaining)*100:.0f}%), tilt={tilt}")

# ============================================================
# Section 7: Statistical tests
# ============================================================

def section7(sites, quick=False):
    print("\n" + "="*60)
    print("SECTION 7: STATISTICAL TESTS")
    print("="*60)
    
    s7 = [s for s in sites if s["combined"] >= 7]
    xyzs_s7 = np.array([ll2xyz(s["lat"], s["lon"]) for s in s7])
    
    tilt_opt, actual_500, pole_opt = optimize_500(s7)
    print(f"\nOptimized pole for score >= 7: tilt={tilt_opt}, captures {actual_500}/{len(s7)}")
    
    print(f"\nCompetitor circles (score >= 7, full scan, {GRID} deg grid):")
    best_counts = {}
    for lat in range(-90, 91, GRID):
        for lon in range(-180, 180, GRID):
            p = ll2xyz(lat, lon)
            dots = np.abs(xyzs_s7 @ p)
            angs = np.arccos(np.clip(dots, 0, 1)) * 180 / np.pi
            best_counts[(lat, lon)] = int(np.sum(np.abs(90 - angs) <= TDEG500))
    
    sorted_poles = sorted(best_counts.items(), key=lambda x: -x[1])
    top_distinct = [sorted_poles[0]]
    for pole, count in sorted_poles[1:]:
        is_distinct = True
        for pp, pc in top_distinct:
            p1 = ll2xyz(pole[0], pole[1])
            p2 = ll2xyz(pp[0], pp[1])
            ang = math.degrees(math.acos(np.clip(np.dot(p1, p2), -1, 1)))
            if ang < 20 or (180 - ang) < 20:
                is_distinct = False
                break
        if is_distinct:
            top_distinct.append((pole, count))
        if len(top_distinct) >= 5:
            break
    
    for pole, count in top_distinct:
        print(f"  ({pole[0]:+d}, {pole[1]:+d}), tilt={90-abs(pole[0])}, count={count}/{len(s7)}")
    
    zero_pct = sum(1 for v in best_counts.values() if v == 0) / len(best_counts) * 100
    print(f"  {zero_pct:.0f}% of poles capture zero sites")
    
    if quick:
        print("\n  [--quick mode: skipping Monte Carlo]")
        return
    
    np.random.seed(42)
    n_trials = 500
    print(f"\nMonte Carlo ({n_trials} trials, latitude-preserving)...")
    
    lats = np.array([s["lat"] for s in s7])
    pole_grid = np.array([ll2xyz(lat, lon) for lat in range(-90, 91, 6) for lon in range(-180, 180, 8)])
    
    mc_counts = []
    for trial in range(n_trials):
        rand_lons = np.random.uniform(-180, 180, len(s7))
        rxyz = np.array([ll2xyz(lats[j], rand_lons[j]) for j in range(len(s7))])
        best = 0
        for p in pole_grid:
            dots = np.abs(rxyz @ p)
            angs = np.arccos(np.clip(dots, 0, 1)) * 180 / np.pi
            c = int(np.sum(np.abs(90 - angs) <= TDEG500))
            if c > best:
                best = c
        mc_counts.append(best)
    
    mc = np.array(mc_counts)
    z = (actual_500 - mc.mean()) / max(mc.std(), 0.01)
    
    print(f"  Actual: {actual_500}")
    print(f"  Null: mean={mc.mean():.2f}, std={mc.std():.2f}, max={mc.max()}")
    print(f"  Z-score: {z:.2f}")
    print(f"  Null >= actual: {np.sum(mc >= actual_500)}/{n_trials}")

# ============================================================
# Main
# ============================================================

def main():
    global GRID
    
    fast = "--fast" in sys.argv
    quick = "--quick" in sys.argv
    section_filter = None
    if "--section" in sys.argv:
        idx = sys.argv.index("--section")
        if idx + 1 < len(sys.argv):
            section_filter = int(sys.argv[idx + 1])
    
    if fast:
        GRID = 4
        print(f"Fast mode: using {GRID} deg grid (results may differ slightly from report)")
    else:
        GRID = 2
        print(f"Using {GRID} deg grid (matches report exactly)")
    
    sites = load_sites()
    print(f"Loaded {len(sites)} sites.")
    
    sections = {
        1: section1, 2: section2, 3: section3,
        4: section4, 5: section5, 6: section6,
    }
    
    if section_filter:
        if section_filter == 7:
            section7(sites, quick=quick)
        elif section_filter in sections:
            sections[section_filter](sites)
    else:
        for num, func in sorted(sections.items()):
            func(sites)
        section7(sites, quick=quick)
    
    print("\n" + "="*60)
    print("COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()
