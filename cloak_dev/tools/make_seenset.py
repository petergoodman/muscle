"""Write the seen-set variant of f3559af's cloak.cpp: test oracle O5 (scripts/test_seenset.sh, tests/README.md).
Usage: make_seenset.py F3559AF_CLOAK_CPP OUT_CPP
Changes, each applied exactly once:
  1. add_to_output_or_discard returns at entry if its argument was already expanded in this MSA-0 column
     (the unused `seen` parameter becomes that set; it is created per column in create_output).
  2. the final quadratic de-duplication becomes an order-preserving hash-set de-duplication.
Both use ndkey(), an injective encoding of a NumOrDash vector under equalNumOrDashVector's equality."""
import sys

src = open(sys.argv[1], encoding='utf-8-sig', newline='').read().replace('\r\n', '\n')
assert '\r' not in src

REPS = [
    ('#include <fstream>\n', '#include <fstream>\n#include <unordered_set>\n'),
    # injective key: length, then (position, type tag, active value) for every entry that is not
    # CHARACTER '-'; an unlisted position is CHARACTER '-'. Equal keys <=> equalNumOrDashVector true.
    ('// Recursively filters and adds valid divvied alignments to the output set.\n',
     'static std::string ndkey(const std::vector<NumOrDash>& v)\n'
     '{\n'
     '    std::string s;\n'
     '    const uint64_t len = v.size();\n'
     '    s.append(reinterpret_cast<const char*>(&len), sizeof(len));\n'
     '    for (size_t i = 0; i < v.size(); ++i) {\n'
     '        const NumOrDash& e = v[i];\n'
     '        if (e.type == NumOrDashType::CHARACTER && e.charValue == \'-\')\n'
     '            continue;\n'
     '        const uint64_t pos = i;\n'
     '        s.append(reinterpret_cast<const char*>(&pos), sizeof(pos));\n'
     '        if (e.type == NumOrDashType::INTEGER) {\n'
     '            s.push_back(\'I\');\n'
     '            s.append(reinterpret_cast<const char*>(&e.intValue), sizeof(int));\n'
     '        } else {\n'
     '            s.push_back(\'C\');\n'
     '            s.push_back(e.charValue);\n'
     '        }\n'
     '    }\n'
     '    return s;\n'
     '}\n'
     '\n'
     '// Recursively filters and adds valid divvied alignments to the output set.\n'),
    ('    const std::vector<NumOrDash>& toSearch,\n'
     '    std::vector<std::vector<NumOrDash>>& seen)\n'
     '{\n',
     '    const std::vector<NumOrDash>& toSearch,\n'
     '    std::unordered_set<std::string>& seen)\n'
     '{\n'
     '    if (!seen.insert(ndkey(toSearch)).second)\n'
     '        return;   // already expanded in this MSA-0 column: would append only repeats\n'),
    ('            std::vector<std::vector<NumOrDash>> seen; // initially empty\n',
     '            std::unordered_set<std::string> seen; // divvies expanded in this column\n'),
    ('    std::vector<std::vector<NumOrDash>> nums_array_output;\n'
     '    for (const auto& num_array : nums_array_output_initial)\n'
     '    {\n'
     '        bool duplicate = false;\n'
     '        for (const auto& existing : nums_array_output)\n'
     '        {\n'
     '            if (equalNumOrDashVector(existing, num_array))\n'
     '            {\n'
     '                duplicate = true;\n'
     '                break;\n'
     '            }\n'
     '        }\n'
     '        if (!duplicate)\n'
     '            nums_array_output.push_back(num_array);\n'
     '    }\n',
     '    std::vector<std::vector<NumOrDash>> nums_array_output;\n'
     '    {\n'
     '        std::unordered_set<std::string> kept;\n'
     '        for (const auto& num_array : nums_array_output_initial)\n'
     '            if (kept.insert(ndkey(num_array)).second)\n'
     '                nums_array_output.push_back(num_array);\n'
     '    }\n'),
]

for old, new in REPS:
    c = src.count(old)
    assert c == 1, (c, old)
    src = src.replace(old, new)
with open(sys.argv[2], 'w', encoding='utf-8-sig', newline='') as o:
    o.write(src)
print('wrote', sys.argv[2])
