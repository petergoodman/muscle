"""Write mutants of src/cloak.cpp for mutation testing (scripts/test_mutants.sh, tests/README.md).

Usage: make_mutants.py ORIG_CLOAK_CPP OUTDIR
Writes OUTDIR/<name>/cloak.cpp for every mutant and OUTDIR/MUTANTS.tsv (name, kind, expected
killer). Kinds: control (unmutated: must survive), equivalent (argued to compute the same output:
must survive), killable (must be killed). Every replacement must apply exactly once, so a mutant
cannot silently stop matching the code. M01-M14 are the 2026-10-30 audit's mutants of
create_output; N01-N15 and X01-X03 were added for cloak-v0.3.0.
"""
import os
import sys

SORT = "        std::sort(keyed.begin(), keyed.end());\n"
KEYLOOP = ("            for (size_t k = 1; k < K; ++k)\n"
           "                key[k - 1] = pos_col[k][r][p];\n")
EMIT_OLD = ("            if (j - i >= minimumColumns) {\n"
            "                std::vector<NumOrDash> column(n, NumOrDash('-'));\n"
            "                for (size_t t = i; t < j; ++t)\n"
            "                    column[keyed[t].second] = first_alignment[keyed[t].second][x];\n"
            "                output.push_back(std::move(column));\n"
            "            }\n"
            "            i = j;\n"
            "        }\n")
GRP_DECL = (SORT +
            "        struct G { std::vector<int> key; size_t row0; size_t size; std::vector<NumOrDash> col; };\n"
            "        std::vector<G> grp;\n")
WRITE = "    write_fasta(outputFilename, firstLabels, finalSeqs);\n"
LABELS0 = "    const auto& firstLabels = allLabels[0];\n"
GAP_NUMS = "                if (isgap(letter)) {\n"


def emit_reordered(sortkey):
    """Collect the groups of one MSA-0 column, sort them by `sortkey`, then emit."""
    return ("            if (j - i >= minimumColumns) {\n"
            "                std::vector<NumOrDash> column(n, NumOrDash('-'));\n"
            "                for (size_t t = i; t < j; ++t)\n"
            "                    column[keyed[t].second] = first_alignment[keyed[t].second][x];\n"
            "                grp.push_back({keyed[i].first, keyed[i].second, j - i, std::move(column)});\n"
            "            }\n"
            "            i = j;\n"
            "        }\n"
            "        std::stable_sort(grp.begin(), grp.end(), [](const G& a, const G& b) { return "
            + sortkey + "; });\n"
            "        for (auto& g : grp) output.push_back(std::move(g.col));\n")


def row_filter(threshold):
    """Drop rows with at most `threshold` residues before writing (cloak-v0.1.0's row deletion)."""
    return ("    std::vector<std::string> keptLabels, keptSeqs;\n"
            "    for (size_t r = 0; r < finalSeqs.size(); ++r)\n"
            "        if (std::count_if(finalSeqs[r].begin(), finalSeqs[r].end(),\n"
            "                [](char c) { return c != '-'; }) > " + str(threshold) + ") {\n"
            "            keptLabels.push_back(firstLabels[r]);\n"
            "            keptSeqs.push_back(finalSeqs[r]);\n"
            "        }\n"
            "    write_fasta(outputFilename, keptLabels, keptSeqs);\n")


# name: (kind, expected killer, [(old, new), ...])
MUTANTS = {
    'M00_unmutated': ('control', 'none', []),
    'M01_reverse_group_order': ('killable', 'H02, H03', [
        (SORT, SORT + "        std::reverse(keyed.begin(), keyed.end());\n")]),
    'M02_lastMSA_first': ('killable', 'H02, H03, coverage', [
        (KEYLOOP, "            for (size_t k = 1; k < K; ++k)\n"
                  "                key[K - 1 - k] = pos_col[k][r][p];\n")]),
    'M03_order_by_first_row': ('killable', 'H02', [
        (SORT, GRP_DECL), (EMIT_OLD, emit_reordered("a.row0 < b.row0"))]),
    'M04_ignore_MSA1': ('killable', 'H01', [
        (KEYLOOP, "            for (size_t k = 2; k < K; ++k)\n"
                  "                key[k - 1] = pos_col[k][r][p];\n")]),
    'M05_ignore_lastMSA': ('killable', 'H02', [
        (KEYLOOP, "            for (size_t k = 1; k + 1 < K; ++k)\n"
                  "                key[k - 1] = pos_col[k][r][p];\n")]),
    'M06_mincol_gt': ('killable', 'H01', [
        ("if (j - i >= minimumColumns) {", "if (j - i > minimumColumns) {")]),
    'M07_filter_ge2': ('equivalent', 'none: remove_singletons filters at -mincol again', [
        ("if (j - i >= minimumColumns) {", "if (j - i >= 2) {")]),
    'M08_no_residue_count_assert': ('equivalent', 'none: the loader guarantees equal counts', [
        ("            asserta(pos_col[k][r].size() == residues0[r]);\n", "")]),
    'M09_no_empty_guard': ('killable', 'E06 (under the sanitizers if not before)', [
        ("    if (n == 0 || first_alignment[0].empty())\n        return output;\n", "")]),
    'M10_no_sort_adjacent_runs': ('killable', 'H02, random', [(SORT, "")]),
    'M11_ignore_MSA8': ('killable', 'Pfam fixtures (K = 16)', [
        (KEYLOOP, "            for (size_t k = 1; k < K; ++k)\n"
                  "                if (k != 8) key[k - 1] = pos_col[k][r][p];\n")]),
    'M12_order_by_MSA1col_then_row': ('killable', 'H02, H03', [
        (SORT, GRP_DECL), (EMIT_OLD, emit_reordered(
            "a.key[0] != b.key[0] ? a.key[0] < b.key[0] : a.row0 < b.row0"))]),
    'M13_order_by_size_desc': ('killable', 'H02', [
        (SORT, GRP_DECL), (EMIT_OLD, emit_reordered("a.size > b.size"))]),
    'M14_order_by_MSA12_then_row': ('killable', 'Pfam fixtures (needs K >= 4)', [
        (SORT, GRP_DECL), (EMIT_OLD, emit_reordered(
            "a.key[0] != b.key[0] ? a.key[0] < b.key[0] : "
            "(a.key.size() > 1 && a.key[1] != b.key[1] ? a.key[1] < b.key[1] : a.row0 < b.row0)"))]),
    'N01_row_filter_gt1': ('killable', 'H04-H07', [(WRITE, row_filter(1))]),
    'N02_row_filter_gt0': ('killable', 'H04, H07', [(WRITE, row_filter(0))]),
    'N03_dash_only_in_nums': ('killable', 'H08', [(GAP_NUMS, "                if (letter == '-') {\n")]),
    'N04_dash_only_in_letters': ('killable', 'H08', [
        ("            if (isgap(ch)) {\n", "            if (ch == '-') {\n")]),
    'N05_nonletters_as_gaps': ('killable', 'H13', [
        (GAP_NUMS, "                if (isgap(letter) || !isalpha((unsigned char) letter)) {\n")]),
    'N06_no_zero_msa_die': ('killable', 'E04', [
        ("    if (E.GetMSACount() == 0)   // a path list naming no MSA loads with only a warning\n"
         "        Die(\"cloak: the ensemble %s holds no MSA\", ensemblePath.c_str());\n", "")]),
    'N07_mincol_returns': ('killable', 'E01', [
        ("    if (minimumColumns < 2)\n        Die(\"cloak: -mincol must be >= 2\");\n",
         "    if (minimumColumns < 2)\n        return;\n")]),
    'N08_bad_alloc_returns': ('killable', 'E16', [
        ("        Die(\"cloak: out of memory\");\n", "        return;\n")]),
    'N09_no_survivor_die_removed': ('killable', 'E15', [
        ("    if (nums_array_output.empty())\n"
         "        Die(\"cloak: no column has %u or more residues\", minimumColumns);\n", "")]),
    'N10_write_unchecked': ('killable', 'E20a-c', [
        ("    ok = (toStdout ? fflush(f) == 0 : fclose(f) == 0) && ok;\n",
         "    (void) (toStdout ? fflush(f) : fclose(f));\n    ok = true;\n")]),
    'N11_no_label_sort': ('killable', 'H12, R1', [
        ("        std::sort(rowData.begin(), rowData.end(),\n"
         "            [](const LabelSeq& a, const LabelSeq& b) {\n"
         "                return a.label < b.label;\n"
         "            });\n", "")]),
    'N12_p_not_p_minus_1': ('killable', 'H08, tagging', [
        ("const size_t p = static_cast<size_t>(v.intValue - 1);",
         "const size_t p = static_cast<size_t>(v.intValue);")]),
    'N13_strip_gt_from_labels': ('killable', 'H16', [
        (LABELS0, "    auto firstLabels = allLabels[0];\n"
                  "    for (auto& l : firstLabels)\n"
                  "        if (!l.empty() && l[0] == '>')\n"
                  "            l.erase(0, 1);\n")]),
    'N14_partial_file_kept': ('killable', 'E20c', [
        ("        if (!toStdout && stat(path.c_str(), &st) == 0 && (st.st_mode & S_IFMT) == S_IFREG)\n"
         "            remove(path.c_str());\n", "        (void) st;\n")]),
    'N15_dash_is_a_file': ('killable', 'H17, E20b', [
        ("    const bool toStdout = (path == \"-\");\n", "    const bool toStdout = false;\n")]),
    'X01_stable_sort': ('equivalent', 'none: keys end in the row, so all are distinct', [
        (SORT, "        std::stable_sort(keyed.begin(), keyed.end());\n")]),
    'X02_msa0_column_in_key': ('equivalent', 'none: every residue keyed together shares column x', [
        ("            std::vector<int> key(K - 1);\n",
         "            std::vector<int> key(K);\n            key[K - 1] = static_cast<int>(x);\n")]),
    'X03_labels_from_last_msa': ('equivalent', 'none: every MSA has the same sorted labels', [
        (LABELS0, "    const auto& firstLabels = allLabels.back();\n")]),
}


def main():
    orig_path, outdir = sys.argv[1], sys.argv[2]
    src = open(orig_path, encoding='utf-8-sig', newline='').read().replace('\r\n', '\n')
    rows = []
    for name, (kind, killer, reps) in MUTANTS.items():
        s = src
        for old, new in reps:
            c = s.count(old)
            assert c == 1, (name, c, old)
            s = s.replace(old, new)
        assert (s != src) == (kind != 'control'), name
        d = os.path.join(outdir, name)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'cloak.cpp'), 'w', encoding='utf-8-sig', newline='') as o:
            o.write(s)
        rows.append(f'{name}\t{kind}\t{killer}\n')
    with open(os.path.join(outdir, 'MUTANTS.tsv'), 'w') as o:
        o.write('name\tkind\texpected_killer\n' + ''.join(rows))
    print(f'{len(rows)} mutants written to {outdir}')


if __name__ == '__main__':
    main()
