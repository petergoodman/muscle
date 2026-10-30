/*
 * CLOAK
 * Originally authored by Chiragdeep Chatur (Masel Lab, University of Arizona)
 * Translated into C++ and integrated into MUSCLE by Peter Goodman (Masel Lab, University of Arizona)
 * Description: CLOAK takes in a set of variant alignments of the same sequences and searches for 
 * consensus among them. 
 * To generate the required input alignments for CLOAK, use Muscle5 options to conduct three 
 * perturbations of its core Hidden Markov Model (HMM) and three perturbations of its core 
 * guide tree, the combination of which produces 16 alternative alignments. 
 * The CLOAK algorithm retains each pair of amino acids in the same column if and only if they are 
 * found together (in the same column) in all 16 alternative alignments. 
 * If multiple subsets of amino acids within a column are found together across all alignments, 
 * but the whole column is not, the subsets are split into multiple columns that are all retained.
 *
 * Usage:       muscle -cloak input_ensemble_file -mincol <integer> -output <output_file_name>
 *
 * Arguments:
 *    - input_ensemble_file : Path to the input MSA file, which can either be an EFA file 
 *                            or a text file with paths to individual MSAs on each line
 *    - -mincol <integer>   : Minimum number of non-gap characters required per column
 *                            for that column to be retained in the output.
 *                            Default value of 2 if not specified
 *    - -output <filename>  : Name of the file where the filtered MSA will be written.
 *
 * Notes:
 *    - This tool is intended to be used as a plug-in with the MUSCLE suite.
 *    - Columns not meeting the minimum occupancy threshold will be discarded.
 *    - Input ensemble must either be EFA file or text file with paths to individual MSAs
 *
 * Example:
 *    muscle -cloak input_aligmments.efa -mincol 5 -output cleaned_alignment.fasta
 *
 */

#include "muscle.h"
#include "myutils.h"
#include "ensemble.h"
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <cassert>
#include <fstream>


//Converts the users input MSA ensemble into a 3-D array
//Takes either an EFA of text file
//3-D array format: MSA->Row->Char
std::vector<std::vector<std::vector<char>>> convertEnsembleTo3D(
    const Ensemble& E,
    std::vector<std::vector<std::string>>& allLabels)
{
    // Build a 3D array: each MSA => 2D array => rows of chars
    std::vector<std::vector<std::vector<char>>> allAlignments;

    uint nMSAs = E.GetMSACount();
    allLabels.resize(nMSAs);

    for (uint i = 0; i < nMSAs; i++)
    {
        const MSA& m = E.GetMSA(i);
        uint seqCount = m.GetSeqCount();
        uint colCount = m.GetColCount();

        // Gather (label, sequence) pairs, then sort by label
        struct LabelSeq
        {
            std::string label;
            std::vector<char> seq;
        };
        std::vector<LabelSeq> rowData;
        rowData.reserve(seqCount);

        for (uint s = 0; s < seqCount; s++)
        {
            LabelSeq ls;
            ls.label = m.GetLabel(s);
            ls.seq.reserve(colCount);

            for (uint c = 0; c < colCount; c++)
            {
                ls.seq.push_back(m.GetChar(s, c));
            }
            rowData.push_back(ls);
        }

        // Sort by label
        std::sort(rowData.begin(), rowData.end(),
            [](const LabelSeq& a, const LabelSeq& b) {
                return a.label < b.label;
            });

        // Build 2D char array
        std::vector<std::vector<char>> twoD;
        twoD.reserve(seqCount);

        std::vector<std::string> labelVec;
        labelVec.reserve(seqCount);

        for (auto& ls : rowData)
        {
            twoD.push_back(ls.seq);
            labelVec.push_back(ls.label);
        }

        allAlignments.push_back(twoD);
        allLabels[i] = labelVec;
    }

    return allAlignments; // [MSA_index][row_index][char_index]
}

// Removes leading and trailing whitespace from a string.
std::string trim(const std::string& s) {
    std::string result = s;
    result.erase(result.begin(),
        std::find_if(result.begin(), result.end(),
            [](unsigned char ch) { return !std::isspace(ch); }));
    result.erase(std::find_if(result.rbegin(), result.rend(),
        [](unsigned char ch) { return !std::isspace(ch); }).base(),
        result.end());
    return result;
}


// Helper class/structure/functions to deal with mixed vectors
// We need to convert sequence vectors from letters to ints
// Dashes must stay as Chars, so we use a helper class
enum class NumOrDashType { INTEGER, CHARACTER };

struct NumOrDash {
    NumOrDashType type;
    union {
        int intValue;
        char charValue;
    };

    // Actual constructor for int
    NumOrDash(int v) {
        type = NumOrDashType::INTEGER;
        intValue = v;
    }

    // Actual constructor for char
    NumOrDash(char c) {
        type = NumOrDashType::CHARACTER;
        charValue = c;
    }
};

//   For each alignment (a 2-D vector of characters) in the 3-D array,
//   and for each sequence (row) in that alignment, assign a running count
//   (starting at 1) to each non-dash character, leaving dashes unchanged.
//   The result is returned as a 3-D vector, where each element is a NumOrDash.
std::vector<std::vector<std::vector<NumOrDash>>>
convert_to_nums(const std::vector<std::vector<std::vector<char>>>& array)
{
    std::vector<std::vector<std::vector<NumOrDash>>> nums;

    // Iterate over each alignment (each 2-D array)
    for (size_t x = 0; x < array.size(); ++x) {
        const auto& two_d_array = array[x];
        std::vector<std::vector<NumOrDash>> make_two_d_array;

        // Process each sequence (row) in the alignment.
        for (size_t y = 0; y < two_d_array.size(); ++y) {
            const auto& row = two_d_array[y];
            int i = 1;  // Reset counter for each sequence.
            std::vector<NumOrDash> row_array;

            // Process each character in the sequence.
            for (size_t z = 0; z < row.size(); ++z) {
                char letter = row[z];
                if (letter == '-') {
                    // Calls the NumOrDash(char) constructor
                    row_array.push_back(NumOrDash('-'));
                }
                else {
                    // Calls the NumOrDash(int) constructor
                    row_array.push_back(NumOrDash(i));
                    i++;
                }
            }
            make_two_d_array.push_back(row_array);
        }
        nums.push_back(make_two_d_array);
    }
    return nums;
}


// Builds the numeric CLOAK output: one entry per retained column, each a vector over rows that
// holds a residue number (1-based, as made by convert_to_nums) or '-'.
// Two residues of a column of alignment 0 are kept together iff they share a column in every
// alignment. Each residue is keyed by the columns it occupies in alignments 1..K-1, so the
// retained columns are exactly the groups of residues with equal keys. Groups are emitted in
// lexicographic key order, which is the order in which the recursive search used up to f3559af
// first reached them, so the output file is unchanged; that search reached each group many times
// over and stored every copy, which exhausted memory on ensembles whose alignments disagree.
// Groups with fewer than minimumColumns residues are dropped here rather than after conversion
// to letters.
std::vector<std::vector<NumOrDash>> create_output(
    const std::vector<std::vector<std::vector<NumOrDash>>>& all_als,
    uint minimumColumns)
{
    std::vector<std::vector<NumOrDash>> output;
    const size_t K = all_als.size();
    const auto& first_alignment = all_als[0];
    const size_t n = first_alignment.size();
    if (n == 0 || first_alignment[0].empty())
        return output;

    // pos_col[k][r][p-1] = column of residue number p of row r in alignment k (k >= 1).
    std::vector<size_t> residues0(n, 0);
    for (size_t r = 0; r < n; ++r)
        for (const auto& v : first_alignment[r])
            if (v.type == NumOrDashType::INTEGER)
                ++residues0[r];
    std::vector<std::vector<std::vector<int>>> pos_col(K);
    for (size_t k = 1; k < K; ++k) {
        pos_col[k].resize(n);
        for (size_t r = 0; r < n; ++r) {
            const auto& row = all_als[k][r];
            for (size_t c = 0; c < row.size(); ++c)
                if (row[c].type == NumOrDashType::INTEGER)
                    pos_col[k][r].push_back(static_cast<int>(c));
            // Equal counts are guaranteed when gaps are '-'; '.' is numbered as a residue here.
            if (pos_col[k][r].size() != residues0[r])
                Die("cloak: row %u has %u residues in MSA %u but %u in MSA 0; gaps must be '-'",
                    (uint) r, (uint) pos_col[k][r].size(), (uint) k, (uint) residues0[r]);
        }
    }

    const size_t num_cols = first_alignment[0].size();
    std::vector<std::pair<std::vector<int>, size_t>> keyed;   // (key, row)
    for (size_t x = 0; x < num_cols; ++x) {
        keyed.clear();
        for (size_t r = 0; r < n; ++r) {
            const NumOrDash& v = first_alignment[r][x];
            if (v.type != NumOrDashType::INTEGER)
                continue;
            const size_t p = static_cast<size_t>(v.intValue - 1);
            std::vector<int> key(K - 1);
            for (size_t k = 1; k < K; ++k)
                key[k - 1] = pos_col[k][r][p];
            keyed.emplace_back(std::move(key), r);
        }
        std::sort(keyed.begin(), keyed.end());
        for (size_t i = 0; i < keyed.size(); ) {
            size_t j = i + 1;
            while (j < keyed.size() && keyed[j].first == keyed[i].first)
                ++j;
            if (j - i >= minimumColumns) {
                std::vector<NumOrDash> column(n, NumOrDash('-'));
                for (size_t t = i; t < j; ++t)
                    column[keyed[t].second] = first_alignment[keyed[t].second][x];
                output.push_back(std::move(column));
            }
            i = j;
        }
    }
    return output;
}


// Helper function, takes in a 2D Array, transposes it
template <typename Container>
std::vector<std::vector<typename Container::value_type>> transpose_output(const std::vector<Container>& l1) {
    std::vector<std::vector<typename Container::value_type>> l2;
    if (l1.empty()) return l2;
    size_t cols = l1[0].size();
    for (size_t i = 0; i < cols; ++i) {
        std::vector<typename Container::value_type> row;
        for (const auto& item : l1) {
            if (i < item.size()) {
                row.push_back(item[i]);
            }
            else {
                break;
            }
        }
        l2.push_back(row);
    }
    return l2;
}


// This function takes a vector of sequences (strings) and returns a new vector 
// containing only those sequences that have more than one non-dash character.
std::vector<std::vector<char>> remove_singletons(
    const std::vector<std::vector<char>>& two_d_array,
    uint minColumns)
{
    std::vector<std::vector<char>> return_array;
    for (const auto& one_d_array : two_d_array) {
        unsigned int num_letters = 0;
        for (char letter : one_d_array) {
            if (letter != '-') {
                ++num_letters;
            }
        }
        if (num_letters >= minColumns) {
            return_array.push_back(one_d_array);
        }
    }
    return return_array;
}


// nums_array: a 2D vector (each inner vector is a row) of NumOrDash (each element is either an int or '-' as a char)
// first_alignment: a 2D vector of char (all letters)
// Returns a 2D vector of char.
std::vector<std::vector<char>> produce_output_in_letters(
    const std::vector<std::vector<NumOrDash>>& nums_array,
    const std::vector<std::vector<char>>& first_alignment)
{
    std::vector<std::vector<char>> return_array;
    std::vector<std::vector<char>> intermediary;

    // If first_alignment is empty, there's nothing to do
    if (first_alignment.empty() || first_alignment[0].empty()) {
        return return_array;
    }

    // Build an "intermediary" 2D array that collects the non-dash
    // characters column by column from first_alignment.
    size_t num_cols = first_alignment[0].size();
    for (size_t i = 0; i < num_cols; i++) {
        std::string s;
        // Concatenate the i-th element from every row
        for (const auto& row : first_alignment) {
            if (i < row.size()) {
                s.push_back(row[i]);
            }
        }
        // Build an array of characters from s, skipping dashes.
        std::vector<char> some_array;
        for (char ch : s) {
            if (ch == '-') {
                continue;
            }
            some_array.push_back(ch);
        }
        intermediary.push_back(some_array);
    }

    // For each row in nums_array:
    for (const auto& nums : nums_array) {
        std::vector<char> subarray;
        size_t cols = nums.size();

        for (size_t i = 0; i < cols; i++) {
            // If it is a dash, we output '-'.
            if (nums[i].type == NumOrDashType::CHARACTER &&
                nums[i].charValue == '-')
            {
                subarray.push_back('-');
            }
            // Else if it's an integer:
            else if (nums[i].type == NumOrDashType::INTEGER) {
                int int_val = nums[i].intValue;
                int index = int_val - 1;  // old code: int index = std::get<int>(nums[i]) - 1;

                // Check bounds in intermediary[i]
                if (i < intermediary.size() && index >= 0 &&
                    index < static_cast<int>(intermediary[i].size()))
                {
                    subarray.push_back(intermediary[i][index]);
                }
                else {
                    // If out-of-range, append a dash.
                    subarray.push_back('-');
                }
            }
            else {
                // Fallback case
                subarray.push_back('-');
            }
        }
        return_array.push_back(subarray);
    }
    return return_array;
}


void cmd_cloak()

{
    // Parse user options
    std::string ensemblePath = opt(cloak);  // user runs: muscle -cloak myfile ...
    if (ensemblePath.empty())
    {
        Log("Error: you must specify -cloak <ensembleFile or list-of-FASTA>\n");
        return;
    }

    std::string outputFilename = opt(output);
    if (outputFilename.empty())
        outputFilename = "cloak_result.fasta";

    uint minimumColumns = 2; // default if -mincol is not specified
    uint userMincol = opt(mincol);  // Parser returns 0 if not specified

    // If the parser's sentinel is 0 when -mincol isn't used, then check:
    if (userMincol != 0)
    {
        // User actually specified -mincol
        minimumColumns = userMincol;
    }
    if (minimumColumns < 2)
    {
        Log("Error: -mincol must be >= 2\n");
        return;
    }

    Log("Starting cloak. Ensemble path: %s, Output filename: %s, Minimum columns requested: %u\n",
        ensemblePath.c_str(), outputFilename.c_str(), userMincol);

    // Load the entire ensemble from the user-supplied file
    Ensemble E;
    try
    {
        E.FromFile(ensemblePath);
    }
    catch (...)
    {
        Log("Error: Unable to parse ensemble from file %s\n", ensemblePath.c_str());
        return;
    }

    // Convert the ensemble into a 3D array of chars
    // plus keep a parallel 2D array of labels for each MSA.
    // allLabels[i][r] = label for row r in the i-th MSA

    Log("Converting ensemble to 3D array of letters...\n");

    std::vector<std::vector<std::string>> allLabels;
    auto three_d_array_letters = convertEnsembleTo3D(E, allLabels);


    // CLOAK Algorithm
    // 1) Convert letters -> numbers (and dashes)
    Log("Converting letters to numeric representation.\n");
    auto nums_answer = convert_to_nums(three_d_array_letters);

    // 2) Create numeric output (each retained column once, so no duplicate removal is needed)
    Log("Creating initial numeric output.\n");
    auto nums_array_output = create_output(nums_answer, minimumColumns);

    // 4) Use the first alignment in three_d_array_letters as reference.
    auto reference_alignment = three_d_array_letters[0];

    // 5) Convert numeric output back to letters
    Log("Converting numeric output back to letters.\n");
    auto transposed_reference = transpose_output(reference_alignment);
    auto letters_array_output_initial =
        produce_output_in_letters(nums_array_output, transposed_reference);

    Log("Removing singletons with minimum column threshold = %u.\n", minimumColumns);
    auto letters_array_output = remove_singletons(letters_array_output_initial, minimumColumns);

    letters_array_output = transpose_output(letters_array_output);

    // Build final MSA with original labels from the original MSAs

    // allLabels[0] = the labels (in alphabetical order) for the first MSA
    const auto& firstLabels = allLabels[0];

    // Build finalLabels/finalSeqs in parallel
    std::vector<std::string> finalLabels;
    std::vector<std::string> finalSeqs;

    for (size_t i = 0; i < letters_array_output.size(); i++)
    {
        // Convert the row to a string
        std::string rowStr(letters_array_output[i].begin(),
            letters_array_output[i].end());
        // Count non-dash letters
        int letterCount = 0;
        for (char c : rowStr)
            if (c != '-') letterCount++;
        // If this row is NOT a singleton, we keep it
        if (letterCount > 1)
        {
            // The label for row i is firstLabels[i]
            if (i < firstLabels.size())
                finalLabels.push_back(firstLabels[i]);
            else
                finalLabels.push_back(std::string(">Seq") + std::to_string(i + 1));

            finalSeqs.push_back(rowStr);
        }
    }

    Log("Building final MSA with %u sequences.\n", (unsigned)finalSeqs.size());

    // Create a new MSA from these labels & sequences
    MSA mCloaked;
    // We must strip the '>' from labels, because MSA::FromStrings2 expects them *without* a '>' prefix.
    for (auto& lab : finalLabels)
    {
        if (!lab.empty() && lab[0] == '>')
            lab.erase(0, 1);
    }
    mCloaked.FromStrings2(finalLabels, finalSeqs);

    mCloaked.ToFASTAFile(outputFilename);

    Log("Cloak: Output written to %s\n", outputFilename.c_str());
}