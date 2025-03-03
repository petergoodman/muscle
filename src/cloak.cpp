#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <variant>
#include <cassert>
#include "myutils.h"
#include "myopts.h"
#include "msa.h"
#include "ensemble.h"

// For convenience.
namespace fs = std::filesystem;

std::vector<std::vector<std::vector<char>>> convertEnsembleTo3D(
    const Ensemble& E,
    std::vector<std::vector<std::string>>& allLabels)
{
    // We'll build a 3D array: each MSA => 2D array => rows of chars
    std::vector<std::vector<std::vector<char>>> allAlignments;

    uint nMSAs = E.GetMSACount();
    allLabels.resize(nMSAs);

    for (uint i = 0; i < nMSAs; i++)
    {
        const MSA& m = E.GetMSA(i);  // i-th MSA
        uint seqCount = m.GetSeqCount();
        uint colCount = m.GetColCount();

        // We'll gather (label, sequence) pairs, then sort by label
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
            ls.label = m.GetLabel(s); // e.g. "species1"
            ls.seq.reserve(colCount);

            for (uint c = 0; c < colCount; c++)
            {
                ls.seq.push_back(m.GetChar(s, c));
            }
            rowData.push_back(ls);
        }

        // Sort rowData by label
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
        allLabels[i] = labelVec; // store the sorted labels
    }

    return allAlignments; // [MSA_index][row_index][char_index]
}

// Helper function: trim
// Removes leading and trailing whitespace from a string.
std::string trim(const std::string& s) {
    std::string result = s;
    // left trim
    result.erase(result.begin(),
        std::find_if(result.begin(), result.end(),
            [](unsigned char ch) { return !std::isspace(ch); }));
    // right trim
    result.erase(std::find_if(result.rbegin(), result.rend(),
        [](unsigned char ch) { return !std::isspace(ch); }).base(),
        result.end());
    return result;
}

// Type alias for an element that can be either an integer (number) or a dash (char).
using NumOrDash = std::variant<int, char>;

// Function: convert_to_nums
// Description:
//   For each alignment (a 2-D vector of characters) in the 3-D array,
//   and for each sequence (row) in that alignment, assign a running count
//   (starting at 1) to each non-dash character, leaving dashes unchanged.
//   The result is returned as a 3-D vector, where each element is a NumOrDash.
std::vector<std::vector<std::vector<NumOrDash>>>
convert_to_nums(const std::vector<std::vector<std::vector<char>>>& array) {
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
                    row_array.push_back('-');  // Leave dashes unchanged.
                }
                else {
                    row_array.push_back(i);    // Replace letter with its number.
                    i++;
                }
            }
            make_two_d_array.push_back(row_array);
        }
        nums.push_back(make_two_d_array);
    }
    return nums;
}

// Function: can_this_divvy_be_found_without_subsequent_divvies
// Parameters:
//   - toSearch: A 1-D vector of NumOrDash representing the column to match.
//   - three_d_array_of_divvies: A 3-D vector, where each element is an alignment
//       (a 2-D vector of rows) and each row is a vector of NumOrDash.
//   - current_alignment_number: The index of the alignment that should be skipped.
// Operation:
//   For each alignment (except the one specified by current_alignment_number),
//   the function checks whether there is at least one row (divvy) such that for every
//   index k, if toSearch[k] is not a dash then toSearch[k] equals divvy[k].
// Returns:
//   true if every alignment (other than the current one) has at least one matching row;
//   otherwise, false.
bool can_this_divvy_be_found_without_subsequent_divvies(
    const std::vector<NumOrDash>& toSearch,
    const std::vector<std::vector<std::vector<NumOrDash>>>& three_d_array_of_divvies,
    size_t current_alignment_number)
{
    // Iterate over each alignment.
    for (size_t x = 0; x < three_d_array_of_divvies.size(); ++x) {
        if (x == current_alignment_number)
            continue;  // Skip the current alignment.
        const auto& divvies_of_current_alignment = three_d_array_of_divvies[x];
        bool current_alignment_has_match = false;
        // For each row (divvy) in this alignment...
        for (const auto& divvy : divvies_of_current_alignment) {
            bool match_was_found = true;
            // Compare each element in the toSearch vector.
            for (size_t k = 0; k < toSearch.size(); ++k) {
                const NumOrDash& current_num = toSearch[k];
                const NumOrDash& compare_to = divvy[k];
                // If current_num is a dash, skip the comparison.
                if (current_num.index() == 1 && std::get<char>(current_num) == '-') {
                    continue;
                }
                else {
                    // Otherwise, first ensure the types match.
                    if (current_num.index() != compare_to.index()) {
                        match_was_found = false;
                        break;
                    }
                    // Compare the values.
                    if (current_num.index() == 0) { // int type
                        if (std::get<int>(current_num) != std::get<int>(compare_to)) {
                            match_was_found = false;
                            break;
                        }
                    }
                    else { // char type
                        if (std::get<char>(current_num) != std::get<char>(compare_to)) {
                            match_was_found = false;
                            break;
                        }
                    }
                }
            } // end for each element k
            if (match_was_found) {
                current_alignment_has_match = true;
                break;
            }
        } // end for each divvy
        if (!current_alignment_has_match) {
            return false;
        }
    }
    return true;
}

// Function: divvies_as_3darray
// For each dictionary in the input array (each representing one alignment),
// iterate over each key (each with an associated vector of NumOrDash).
// Count the non-dash elements in that vector, and if there are at least 2,
// include that vector in the alignment's output. Finally, return a 3-D array
// (vector of alignments, where each alignment is a vector of rows, and each row
// is a vector of NumOrDash).                  
std::vector<std::vector<std::vector<NumOrDash>>>
divvies_as_3darray(const std::vector<std::map<int, std::vector<NumOrDash>>>& array_of_dicts) {
    std::vector<std::vector<std::vector<NumOrDash>>> return_array;

    // Iterate over each dictionary (each alignment)
    for (const auto& divvy_dict : array_of_dicts) {
        std::vector<std::vector<NumOrDash>> array_per_alignment;

        // For each key-value pair in the dictionary:
        for (const auto& entry : divvy_dict) {
            const auto& curr_array = entry.second;
            int number_of_nums = 0;
            // Count the number of elements that are not dashes.
            for (const auto& elem : curr_array) {
                // If elem holds an int, it is automatically not a dash.
                if (elem.index() == 0) {
                    number_of_nums++;
                }
                else {
                    // If it holds a char, check if it is not '-'.
                    if (std::get<char>(elem) != '-') {
                        number_of_nums++;
                    }
                }
            }
            // If at least 2 non-dash elements exist, include this vector.
            if (number_of_nums >= 2) {
                array_per_alignment.push_back(curr_array);
            }
        }
        return_array.push_back(array_per_alignment);
    }

    return return_array;
}

// Function: get_divvied_results_as_array_of_dicts
// For each indices array in all_indices, build a dictionary (map) where:
//   - Each key is a non-dash number found in the indices array.
//   - Then for each position in the indices array:
//       * If the element is a dash, append '-' to every key’s vector.
//       * Otherwise, let curr be the number; then append to myDict[curr] the
//         corresponding element from toSearch, and for every other key append '-'.
// Returns a vector of such dictionaries.
std::vector<std::map<int, std::vector<NumOrDash>>>
get_divvied_results_as_array_of_dicts(const std::vector<std::vector<NumOrDash>>& all_indices,
    const std::vector<NumOrDash>& toSearch) {
    std::vector<std::map<int, std::vector<NumOrDash>>> someArray;

    // Process each indices array.
    for (const auto& indices : all_indices) {
        std::map<int, std::vector<NumOrDash>> myDict;
        // First pass: add a key for every non-dash element.
        for (size_t i = 0; i < indices.size(); i++) {
            // If indices[i] is not a dash, it should hold an int.
            if (indices[i].index() == 0) { // index 0 corresponds to int
                int key = std::get<int>(indices[i]);
                if (myDict.find(key) == myDict.end()) {
                    myDict[key] = std::vector<NumOrDash>();  // initialize an empty vector.
                }
            }
        }
        // Second pass: fill in the vectors.
        for (size_t i = 0; i < indices.size(); i++) {
            // If the element is a dash, append '-' to every key's vector.
            if (indices[i].index() == 1 && std::get<char>(indices[i]) == '-') {
                for (auto& pair : myDict) {
                    pair.second.push_back(NumOrDash{ '-' });
                }
            }
            else {
                // Otherwise, the element is an int.
                int curr = std::get<int>(indices[i]);
                // Append the corresponding toSearch element for the matching key.
                myDict[curr].push_back(toSearch[i]);
                // For every other key, append '-'.
                for (auto& pair : myDict) {
                    if (pair.first != curr) {
                        pair.second.push_back(NumOrDash{ '-' });
                    }
                }
            }
        }
        someArray.push_back(myDict);
    }
    return someArray;
}

// Function: divvy_is_required
// Description:
//   For each indices array in 'all_indices', this function skips any leading
//   dashes, then takes the first non-dash value as a reference number. It then
//   examines the remainder of the array (ignoring dashes). If any non-dash value
//   is different from the reference, it returns true (divvy is required).
//   If no indices array shows a discrepancy, it returns false.
bool divvy_is_required(const std::vector<std::vector<NumOrDash>>& all_indices) {
    for (const auto& indices_array : all_indices) {
        if (indices_array.empty())
            continue;
        size_t i = 0;
        // Skip over any leading dashes.
        while (i < indices_array.size() &&
            std::holds_alternative<char>(indices_array[i]) &&
            std::get<char>(indices_array[i]) == '-') {
            i++;
        }
        // If the entire array is dashes, skip to the next indices array.
        if (i == indices_array.size()) {
            continue;
        }
        // The first non-dash element should be an int.
        int number_to_match_with = 0;
        if (std::holds_alternative<int>(indices_array[i])) {
            number_to_match_with = std::get<int>(indices_array[i]);
        }
        else {
            // This branch should not occur if our data is as expected.
            continue;
        }
        // Check the remainder of the indices_array.
        for (size_t x = i; x < indices_array.size(); ++x) {
            // If the element is a dash, skip it.
            if (std::holds_alternative<char>(indices_array[x]) &&
                std::get<char>(indices_array[x]) == '-') {
                continue;
            }
            else if (std::holds_alternative<int>(indices_array[x])) {
                int currentNum = std::get<int>(indices_array[x]);
                if (currentNum != number_to_match_with) {
                    return true; // A discrepancy is found.
                }
            }
        }
    }
    return false;
}

// Function: create_indices
// Description:
//   Given the 3-D array of alignments (all_alignments), a 1-D array (toSearch)
//   containing one element from each sequence of the first alignment, and an integer
//   (curr_alignment_number) indicating which alignment to process, this function
//   searches each row of the specified alignment for the element from toSearch.
//   If the element is a dash ('-'), a dash is appended to the result;
//   otherwise, for each occurrence in the row equal to the element, the index is appended.
std::vector<NumOrDash> create_indices(
    const std::vector<std::vector<std::vector<NumOrDash>>>& all_alignments,
    const std::vector<NumOrDash>& toSearch,
    size_t curr_alignment_number)
{
    std::vector<NumOrDash> indices;
    const auto& current_alignment = all_alignments[curr_alignment_number];

    for (size_t x = 0; x < current_alignment.size(); ++x) {
        const auto& row = current_alignment[x];
        NumOrDash num = toSearch[x];

        // If the toSearch element is a dash, record a dash.
        if (std::holds_alternative<char>(num) && std::get<char>(num) == '-') {
            indices.push_back('-');
        }
        else {
            // Otherwise, for each element in the row, if it equals num, append its index.
            for (size_t a = 0; a < row.size(); ++a) {
                // Compare only if the types match.
                if (row[a].index() == num.index()) {
                    if (row[a].index() == 0) {  // int type
                        if (std::get<int>(row[a]) == std::get<int>(num)) {
                            indices.push_back(static_cast<int>(a));
                        }
                    }
                    else {  // char type
                        if (std::get<char>(row[a]) == std::get<char>(num)) {
                            indices.push_back(static_cast<int>(a));
                        }
                    }
                }
            }
        }
    }
    return indices;
}

std::vector<NumOrDash> create_toSearch(const std::vector<std::vector<NumOrDash>>& alignment, size_t col) {
    std::vector<NumOrDash> result;
    for (const auto& row : alignment) {
        if (col < row.size())
            result.push_back(row[col]);
    }
    return result;
}


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
std::vector<std::string> remove_singletons(const std::vector<std::string>& two_d_array) {
    std::vector<std::string> return_array;
    for (const std::string& one_d_array : two_d_array) {
        int num_letters = 0;
        for (char letter : one_d_array) {
            if (letter != '-') {
                ++num_letters;
            }
        }
        if (num_letters > 1) {
            return_array.push_back(one_d_array);
        }
    }
    return return_array;
}

void add_to_output_or_discard(
    const std::vector<std::vector<NumOrDash>>& indices_of_all_alignments,
    std::vector<std::vector<NumOrDash>>& output,
    const std::vector<NumOrDash>& toSearch,
    std::vector<std::vector<NumOrDash>>& seen)
{
    // Count the number of non-dash elements in toSearch.
    int number_of_nums = 0;
    for (const auto& num : toSearch) {
        if (num.index() == 0) { // int
            number_of_nums++;
        }
        else if (std::get<char>(num) != '-') {
            number_of_nums++;
        }
    }
    if (number_of_nums == 0) {
        return;
    }

    auto array_containing_dicts_of_divvies = get_divvied_results_as_array_of_dicts(indices_of_all_alignments, toSearch);
    auto three_d_array_of_divvies = divvies_as_3darray(array_containing_dicts_of_divvies);

    // Loop starting from index 1 (skip alignment 0).
    for (size_t i = 1; i < three_d_array_of_divvies.size(); i++) {
        const auto& divvies_from_specific_alignment = three_d_array_of_divvies[i];
        for (const auto& divvy : divvies_from_specific_alignment) {
            // If toSearch equals divvy, skip it.
            if (divvy == toSearch)
                continue;
            bool can_found = can_this_divvy_be_found_without_subsequent_divvies(divvy, three_d_array_of_divvies, i);
            if (can_found) {
                output.push_back(divvy);
                continue;
            }
            else {
                add_to_output_or_discard(indices_of_all_alignments, output, divvy, seen);
            }
        }
    }
}

//---------------------------------------------------------------------
// Main function: create_output
//---------------------------------------------------------------------
std::vector<std::vector<NumOrDash>> create_output(
    const std::vector<std::vector<std::vector<NumOrDash>>>& all_als)
{
    std::vector<std::vector<NumOrDash>> output;
    const auto& first_alignment = all_als[0];
    // Use the number of elements in the first row as the number of columns.
    size_t num_cols = first_alignment[0].size();

    for (size_t x = 0; x < num_cols; x++) {
        // Create the initial toSearch term (i.e. column x of first_alignment)
        std::vector<NumOrDash> toSearch = create_toSearch(first_alignment, x);

        // Build indices_of_all_alignments: one indices array per alignment.
        std::vector<std::vector<NumOrDash>> indices_of_all_alignments;
        for (size_t k = 0; k < all_als.size(); k++) {
            std::vector<NumOrDash> indices = create_indices(all_als, toSearch, k);
            indices_of_all_alignments.push_back(indices);
        }

        // Check if divvy is needed.
        bool divvy_needed = divvy_is_required(indices_of_all_alignments);
        if (!divvy_needed) {
            output.push_back(toSearch);
        }
        else {
            std::vector<std::vector<NumOrDash>> seen; // initially empty
            add_to_output_or_discard(indices_of_all_alignments, output, toSearch, seen);
        }
        // Reset toSearch is not necessary.
    }
    return output;
}

// produce_output_in_letters:
//  - output: a 2D vector (each inner vector is a row) of NumOrDash (each element is either an int or '-' as a char)
//  - first_alignment: a 2D vector of char (all letters)
// Returns a 2D vector of char.
std::vector<std::vector<char>> produce_output_in_letters(
    const std::vector<std::vector<NumOrDash>>& output,
    const std::vector<std::vector<char>>& first_alignment)
{
    std::vector<std::vector<char>> return_array;
    std::vector<std::vector<char>> intermediary;

    // For each column index in first_alignment:
    if (first_alignment.empty() || first_alignment[0].empty()) {
        return return_array;
    }
    size_t num_cols = first_alignment[0].size();
    for (size_t i = 0; i < num_cols; i++) {
        std::string s = "";
        // Concatenate the i-th element from every row
        for (const auto& row : first_alignment) {
            if (i < row.size())
                s.push_back(row[i]);
        }
        // Build an array of characters from s, skipping dashes.
        std::vector<char> some_array;
        for (char ch : s) {
            if (ch == '-') continue;
            some_array.push_back(ch);
        }
        intermediary.push_back(some_array);
    }

    // For each row in output:
    for (const auto& nums : output) {
        std::vector<char> subarray;
        // Assume every row in output has the same length as intermediary.size()
        size_t cols = nums.size();
        for (size_t i = 0; i < cols; i++) {
            // curr_num is the element from the output row at column i.
            // If it is a dash, we output '-'.
            if (nums[i].index() == 1 && std::get<char>(nums[i]) == '-') {
                subarray.push_back('-');
            }
            else if (nums[i].index() == 0) {
                int int_val = std::get<int>(nums[i]);
                int index = int_val - 1; // because Python subtracts 1
                // Check bounds before using intermediary[i]
                if (i < intermediary.size() && index < static_cast<int>(intermediary[i].size())) {
                    subarray.push_back(intermediary[i][index]);
                }
                else {
                    // If out-of-range, append a dash.
                    subarray.push_back('-');
                }
            }
            else {
                // Fallback case.
                subarray.push_back('-');
            }
        }
        return_array.push_back(subarray);
    }
    return return_array;
}

bool equalNumOrDashVector(const std::vector<NumOrDash>& a, const std::vector<NumOrDash>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i].index() != b[i].index()) return false;
        if (a[i].index() == 0) {
            if (std::get<int>(a[i]) != std::get<int>(b[i]))
                return false;
        }
        else {
            if (std::get<char>(a[i]) != std::get<char>(b[i]))
                return false;
        }
    }
    return true;
}

void cmd_cloak()
{
    // (A) === Parse user options ===
    std::string ensemblePath = opt(cloak);  // user runs: muscle -cloak myfile ...
    if (ensemblePath.empty())
    {
        Log("Error: you must specify -cloak <ensembleFile or list-of-FASTA>\n");
        return;
    }

    std::string outputFilename = opt(output);
    if (outputFilename.empty())
        outputFilename = "cloak_result.fasta";

    // (B) === Load the entire ensemble from the user-supplied file ===
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

    // (C) === Convert the ensemble into a 3D array of chars ===
    //     plus keep a parallel 2D array of labels for each MSA.
    //     allLabels[i][r] = label for row r in the i-th MSA
    std::vector<std::vector<std::string>> allLabels;
    auto three_d_array_letters = convertEnsembleTo3D(E, allLabels);

    // (D) === Cloak Steps ===
    // 1) Convert letters -> numbers (and dashes)
    auto nums_answer = convert_to_nums(three_d_array_letters);

    // 2) Create initial numeric output
    auto nums_array_output_initial = create_output(nums_answer);

    // 3) Remove duplicates
    std::vector<std::vector<NumOrDash>> nums_array_output;
    for (const auto& num_array : nums_array_output_initial)
    {
        bool duplicate = false;
        for (const auto& existing : nums_array_output)
        {
            if (equalNumOrDashVector(existing, num_array))
            {
                duplicate = true;
                break;
            }
        }
        if (!duplicate)
            nums_array_output.push_back(num_array);
    }

    // 4) Use the first alignment in three_d_array_letters as reference (as your old code does).
    //    So reference_alignment = three_d_array_letters[0].
    auto reference_alignment = three_d_array_letters[0];

    // 5) Convert numeric output back to letters
    auto transposed_reference = transpose_output(reference_alignment);
    auto letters_array_output_initial =
        produce_output_in_letters(nums_array_output, transposed_reference);

    // 6) Remove singletons (rows with <=1 letter)
    std::vector<std::string> sequences;
    sequences.reserve(letters_array_output_initial.size());
    for (const auto& row : letters_array_output_initial)
        sequences.push_back(std::string(row.begin(), row.end()));
    sequences = remove_singletons(sequences);

    // 7) Rebuild 2D char array from the pruned sequences
    std::vector<std::vector<char>> letters_array_output;
    letters_array_output.reserve(sequences.size());
    for (auto& s : sequences)
        letters_array_output.push_back(std::vector<char>(s.begin(), s.end()));

    // 8) Transpose again
    letters_array_output = transpose_output(letters_array_output);

    // (E) === Build final MSA with original labels from the *first* MSA ===
    //     Because your code references the "first alignment" as a reference,
    //     we assume the final alignment has the same number of rows as that first MSA
    //     minus any singletons that got removed.

    // allLabels[0] = the labels (in alphabetical order) for the first MSA
    const auto& firstLabels = allLabels[0];

    // The final row count is `letters_array_output.size()`
    // but if the cloak code hasn't changed row ordering, then row i corresponds to label i.
    // We must skip rows that "remove_singletons()" removed.
    // The simplest way is to do the "singleton removal" in parallel with the labels.

    // We'll build finalLabels/finalSeqs in parallel
    std::vector<std::string> finalLabels;
    std::vector<std::string> finalSeqs;

    // letters_array_output is a 2D array: row i => letters_array_output[i]
    // but that "row i" must match firstLabels[i], if it wasn't removed.
    // So let's do something similar to how we removed singletons, but in one step:
    for (size_t i = 0; i < letters_array_output_initial.size(); i++)
    {
        // Convert the row to a string
        std::string rowStr(letters_array_output_initial[i].begin(),
            letters_array_output_initial[i].end());
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

    // Now we have finalLabels and finalSeqs for the final MSA
    // If your code requires that finalSeqs is re-transposed, do that. 
    // But typically we want finalSeqs as row strings, which is what MSA::FromStrings2 expects.
    // At this point, finalSeqs[i] is the i-th row (FASTA row).

    // (F) === Create a new MSA from these labels & sequences ===
    MSA mCloaked;
    // We must strip the '>' from labels, because MSA::FromStrings2 expects them *without* a '>' prefix,
    // typically. Or if your code handles the full string as-is, that might be okay. 
    // We'll assume we remove leading '>' if it exists:
    for (auto& lab : finalLabels)
    {
        if (!lab.empty() && lab[0] == '>')
            lab.erase(0, 1); // remove leading '>'
    }
    mCloaked.FromStrings2(finalLabels, finalSeqs);

    mCloaked.ToFASTAFile(outputFilename);

    Log("Cloak: Output written to %s\n", outputFilename.c_str());
}