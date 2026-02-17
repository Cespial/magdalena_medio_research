#!/usr/bin/env python3
"""
Fix ecosystem_services_results.json to be consistent with final classification areas.
T3 carbon was from OLD classification; must recalculate using NEW T3 areas.
"""
import json
import os

BASE = os.path.dirname(os.path.abspath(__file__))
STATS = os.path.join(BASE, 'outputs', 'phase3_stats')

# IPCC Tier 1 carbon pools (Mg C/ha total)
CARBON = {1: 242, 2: 146, 3: 48.5, 4: 45.5, 5: 0, 6: 22, 7: 15}

# Load data
with open(os.path.join(STATS, 'classification_metrics.json')) as f:
    clf = json.load(f)
with open(os.path.join(STATS, 'ecosystem_services_results.json')) as f:
    eco = json.load(f)

periods = ['pre_acuerdo', 'transicion', 'post_acuerdo_1', 'post_acuerdo_2']
years = [2013, 2016, 2020, 2024]

print("=== Carbon Stock Recalculation ===\n")
print("Using classification_metrics.json areas (FINAL classification):\n")

carbon_stocks = {}
for pk, yr in zip(periods, years):
    areas = clf[pk]['class_areas_ha']
    c_total = 0
    for cid_str, info in areas.items():
        cid = int(cid_str)
        c_total += info['area_ha'] * CARBON.get(cid, 0)
    carbon_stocks[pk] = round(c_total)

    old_c = eco[pk]['carbon_Mg_C']
    diff = c_total - old_c
    print(f"  {pk} ({yr}):")
    print(f"    OLD carbon: {old_c/1e6:.1f} Mt C")
    print(f"    NEW carbon: {c_total/1e6:.1f} Mt C")
    print(f"    Difference: {diff/1e6:+.1f} Mt C")
    print()

# Recalculate carbon changes
changes = {}
pairs = [('pre_acuerdo', 'transicion', '2013_2016'),
         ('transicion', 'post_acuerdo_1', '2016_2020'),
         ('post_acuerdo_1', 'post_acuerdo_2', '2020_2024')]

print("=== Carbon Changes ===\n")
for pk1, pk2, label in pairs:
    change = carbon_stocks[pk2] - carbon_stocks[pk1]
    old_change = eco.get(f'carbon_change_{label}', {}).get('net_Mg_C', 0)
    changes[label] = change
    print(f"  {label}: {change/1e6:+.1f} Mt C (was {old_change/1e6:+.1f} Mt C)")

total_loss = carbon_stocks['post_acuerdo_2'] - carbon_stocks['pre_acuerdo']
print(f"\n  Total 2013-2024: {total_loss/1e6:+.1f} Mt C")
print(f"  Percent change: {total_loss/carbon_stocks['pre_acuerdo']*100:.1f}%")

# Update ecosystem JSON
for pk in periods:
    eco[pk]['carbon_Mg_C'] = float(carbon_stocks[pk])

for pk1, pk2, label in pairs:
    eco[f'carbon_change_{label}']['net_Mg_C'] = float(changes[label])

# Save
out_path = os.path.join(STATS, 'ecosystem_services_results.json')
with open(out_path, 'w') as f:
    json.dump(eco, f, indent=2)
print(f"\n>> Updated: {out_path}")

# Verify
print("\n=== Verification ===")
for pk1, pk2, label in pairs:
    stock_diff = eco[pk2]['carbon_Mg_C'] - eco[pk1]['carbon_Mg_C']
    change = eco[f'carbon_change_{label}']['net_Mg_C']
    match = abs(stock_diff - change) < 1
    print(f"  {label}: stock diff = {stock_diff/1e6:.1f}, change = {change/1e6:.1f}, match = {match}")
