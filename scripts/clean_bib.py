#!/usr/bin/env python3
"""
Clean references.bib to keep only the 54 entries actually cited in the manuscript.
"""

import re
import shutil
from pathlib import Path

# --- Configuration ---
BIB_PATH = Path("/Users/cristianespinal/Documents/magdalena_medio_research/overleaf/references.bib")
CLEAN_PATH = Path("/Users/cristianespinal/Documents/magdalena_medio_research/overleaf/references_clean.bib")

KEEP_KEYS = {
    "Abatzoglou2018",
    "Alvarez2003",
    "Alvarez2012",
    "Anselin1995",
    "Breiman2001",
    "BrizuelaTorres2025Thirty",
    "Cantillo2022Armed",
    "CastroNunez2020Reducing",
    "CastroNunez2022",
    "Chazdon2016",
    "Clerici2020",
    "Congalton1991",
    "Davalos2021Forests",
    "Drusch2012",
    "Fagan2020",
    "Foody2002",
    "Foody2020",
    "Fotheringham2002",
    "Funk2015",
    "Getis1992",
    "Gorelick2017",
    "GutierrezZapata2025Deforestation",
    "Halmy2015",
    "Hansen2013",
    "Hoffmann2018Local",
    "Holdridge1967",
    "IDEAM2018",
    "IPCC2006",
    "Krause2020Reducing",
    "Landholm2019Diverging",
    "Lopez2024Landscape",
    "Moran1950",
    "MurilloSandoval2021",
    "MurilloSandoval2023Disentangling",
    "MurilloSandoval2023Post",
    "Negret2019",
    "Olofsson2014",
    "Oshan2019",
    "Pacheco2009",
    "PerezMarulanda2025Boosting",
    "Pontius2011",
    "Prem2020",
    "Puyravaud2003",
    "Restrepo2006",
    "Rodriguez2023Deforestation",
    "SanchezCuervo2013",
    "Sharp2020",
    "Sierra2017",
    "Stehman2019",
    "Toro2022Interacting",
    "Turner2007",
    "VanegasCubillos2022Forest",
    "Verburg2015",
    "Wulder2022",
}

assert len(KEEP_KEYS) == 54, f"Expected 54 keys, got {len(KEEP_KEYS)}"


def parse_bib_entries(text: str):
    """
    Parse a .bib file into individual entries.
    Handles nested braces correctly.
    Returns a list of (citation_key, full_entry_text) tuples.
    """
    entries = []
    # Pattern to find the start of an entry: @type{key,
    entry_start_pattern = re.compile(r'@(\w+)\s*\{', re.IGNORECASE)

    pos = 0
    while pos < len(text):
        match = entry_start_pattern.search(text, pos)
        if match is None:
            break

        entry_type = match.group(1).lower()
        start = match.start()
        brace_start = match.end() - 1  # position of the opening {

        # Skip @comment, @string, @preamble (they don't have citation keys in the same way)
        if entry_type in ('comment', 'string', 'preamble'):
            # Find the matching closing brace
            depth = 1
            i = brace_start + 1
            while i < len(text) and depth > 0:
                if text[i] == '{':
                    depth += 1
                elif text[i] == '}':
                    depth -= 1
                i += 1
            pos = i
            continue

        # Extract the citation key (everything from after { to the first comma)
        key_start = brace_start + 1
        comma_pos = text.find(',', key_start)
        if comma_pos == -1:
            pos = key_start
            continue

        citation_key = text[key_start:comma_pos].strip()

        # Now find the matching closing brace, counting nesting
        depth = 1
        i = brace_start + 1
        while i < len(text) and depth > 0:
            if text[i] == '{':
                depth += 1
            elif text[i] == '}':
                depth -= 1
            i += 1

        full_entry = text[start:i]
        entries.append((citation_key, full_entry))
        pos = i

    return entries


def main():
    # Read the original file
    raw = BIB_PATH.read_text(encoding='utf-8')
    print(f"Read {BIB_PATH}")
    print(f"File size: {len(raw):,} characters")

    # Parse entries
    entries = parse_bib_entries(raw)
    total = len(entries)
    print(f"Total entries parsed: {total}")

    # Filter
    kept = []
    removed = []
    found_keys = set()

    for key, entry_text in entries:
        if key in KEEP_KEYS:
            kept.append((key, entry_text))
            found_keys.add(key)
        else:
            removed.append(key)

    # Check for missing keys
    missing_keys = KEEP_KEYS - found_keys

    # Build clean output
    clean_text = "\n\n".join(entry_text for _, entry_text in kept) + "\n"

    # Write the clean copy
    CLEAN_PATH.write_text(clean_text, encoding='utf-8')
    print(f"\nWrote clean file to: {CLEAN_PATH}")

    # Overwrite the original
    BIB_PATH.write_text(clean_text, encoding='utf-8')
    print(f"Overwrote original file: {BIB_PATH}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Entries in original file : {total}")
    print(f"  Entries KEPT             : {len(kept)}")
    print(f"  Entries REMOVED          : {len(removed)}")
    print(f"  Target keys provided     : {len(KEEP_KEYS)}")
    print(f"  Target keys found        : {len(found_keys)}")

    if missing_keys:
        print(f"\n  WARNING: {len(missing_keys)} citation key(s) NOT FOUND in the bib file:")
        for k in sorted(missing_keys):
            print(f"    - {k}")
    else:
        print(f"\n  All {len(KEEP_KEYS)} citation keys were found.")

    # List kept keys for verification
    print(f"\nKept entries (alphabetical):")
    for k in sorted(found_keys):
        print(f"  {k}")


if __name__ == "__main__":
    main()
