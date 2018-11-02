// ==========================================================================
//                              SRA_Search - Prototype
// ==========================================================================
// Copyright (c) 2018, Enrico Seiler, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Enrico Seiler or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SEILER OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================

#include <fstream>
#include <unordered_set>
#include <unordered_map>

#include <seqan/arg_parse.h>
#include <seqan/binning_directory.h>

#include "helper.h"

using namespace seqan;

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

struct Options
{
    CharString  contigs_dir;
    CharString  output_file;

    uint32_t    kmer_size;
    uint32_t    window_size;
    uint32_t    number_of_bins;
    unsigned    threads;

    Options():
        kmer_size(19),
        window_size(23),
        number_of_bins(64),
        threads(1) {}
};

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "SRA_search count prototype");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE FILE DIR"));
    setHelpText(parser, 0, "A directory containing reference genome files.");

    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-file", "Specify an output for the counts. \
                                     Default: use the directory name of reference genomes.", ArgParseOption::OUTPUT_FILE));

    addOption(parser, ArgParseOption("b", "number-of-bins", "The number of bins",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "number-of-bins", "1");
    setMaxValue(parser, "number-of-bins", "4194300");

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threads);

    addOption(parser, ArgParseOption("k", "kmer-size", "The size of kmers to count",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "kmer-size", "14");
    setMaxValue(parser, "kmer-size", "32");

    addOption(parser, ArgParseOption("w", "window-size", "The size of the window to count",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "window-size", "14");
}

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse contigs input file.
    getArgumentValue(options.contigs_dir, parser, 0);

    // Append trailing slash if it doesn't exist.
    append_trailing_slash(options.contigs_dir);

    // Parse contigs index prefix.
    getOptionValue(options.output_file, parser, "output-file");
    if (!isSet(parser, "output-file"))
    {
        options.output_file = trimExtension(options.contigs_dir);
        append(options.output_file, "kmer.counts");
    }

    if (isSet(parser, "number-of-bins")) getOptionValue(options.number_of_bins, parser, "number-of-bins");
    if (isSet(parser, "kmer-size")) getOptionValue(options.kmer_size, parser, "kmer-size");
    if (isSet(parser, "window-size")) getOptionValue(options.window_size, parser, "window-size");
    if (isSet(parser, "threads")) getOptionValue(options.threads, parser, "threads");

    return ArgumentParser::PARSE_OK;
}

inline std::unordered_map<uint32_t, std::vector<uint32_t>> map_bins(Options const & options, std::string const & com_ext)
{
    uint32_t batch_size = options.number_of_bins/options.threads;
    if(batch_size * options.threads < options.number_of_bins) ++batch_size;

    std::unordered_map<uint32_t, std::vector<uint32_t>> bin_map;

    std::vector<uint64_t> bin_sizes(options.number_of_bins);

    for (uint32_t bin_number = 0; bin_number < options.number_of_bins; ++bin_number)
    {
        CharString seq_file_path;
        append_file_name(seq_file_path, options.contigs_dir, bin_number);
        append(seq_file_path, com_ext);

        bin_sizes[bin_number] = filesize(toCString(seq_file_path));
    }
    std::cerr << "BIN SIZES\n";
    for (size_t i = 0; i < bin_sizes.size(); ++i)
    {
        std::cerr << i << '\t' << bin_sizes[i];
    }
    std::cerr << '\n';
    uint64_t total_size = std::accumulate(bin_sizes.begin(), bin_sizes.end(), 0);
    double alpha = 1.15;
    double threshold = (1.0 / options.threads) * alpha;
    std::cerr << "Threshold\t" << threshold;
    std::vector<double> bin_weights;
    bin_weights.resize(options.number_of_bins, 0.0);
    std::transform(bin_sizes.begin(), bin_sizes.end(), bin_weights.begin(),
                  [&](auto & size) { return (double) size / total_size;});
    std::cerr << "BIN WEIGHTS\n";
    for (size_t i = 0; i < bin_weights.size(); ++i)
    {
        std::cerr << i << '\t' << bin_weights[i];
    }
    std::cerr << '\n';
    double s{0.0};
    unsigned t{0};
    for (uint32_t bin = 0; bin < options.number_of_bins; ++bin)
    {
        double size_ratio = bin_weights[bin];
        if (s + size_ratio >= threshold)
        {
            s = size_ratio;
            if (t != options.threads -1)
                ++t;
            bin_map[t].push_back(bin);
        }
        else
        {
            s += size_ratio;
            bin_map[t].push_back(bin);
        }
    }
    return bin_map;
}


inline void count_kmers(Options & options)
{
    std::string com_ext = common_ext(options.contigs_dir, options.number_of_bins);

    uint32_t batch_size = options.number_of_bins/options.threads;
    if(batch_size * options.threads < options.number_of_bins) ++batch_size;

    std::unordered_map<uint32_t, std::vector<uint32_t>> bin_map = map_bins(options, com_ext);

    std::vector<std::future<void>> tasks;

    std::mutex print_mtx;

    uint64_t bv_size = 1ULL<<(2*options.kmer_size);
    uint64_t sig_bit = bv_size - 1;

    sdsl::bit_vector overall_content(bv_size, 0, 1);

    for (uint32_t task_number = 0; task_number < options.threads; ++task_number)
    {
        tasks.emplace_back(std::async([=, &print_mtx, &overall_content, &bin_map] {
            for (uint32_t bin_number : bin_map[task_number])
            {
                CharString seq_file_path;
                append_file_name(seq_file_path, options.contigs_dir, bin_number);
                append(seq_file_path, com_ext);

                // read everything as CharString to avoid impure sequences crashing the program
                CharString id, seq;
                SeqFileIn seq_file_in;
                if (!open(seq_file_in, toCString(seq_file_path)))
                {
                    CharString msg = "Unable to open contigs file: ";
                    append(msg, CharString(seq_file_path));
                    std::cerr << msg << std::endl;
                    throw toCString(msg);
                }
                std::unordered_set<uint64_t> hashes;
                while(!atEnd(seq_file_in))
                {
                    BDHash<Dna, Minimizer<19,25>> minimizer;
                    minimizer.resize(options.kmer_size, options.window_size);
                    readRecord(id, seq, seq_file_in);
                    if(length(seq) < options.kmer_size)
                        continue;
                    auto mins = minimizer.getHash(seq);
                    hashes.insert(mins.begin(), mins.end());
                }
                print_mtx.lock();
                std::cerr << bin_number << '\t' << hashes.size() << std::endl;
                print_mtx.unlock();
                for (auto &x : hashes)
                {
                    overall_content[(x & sig_bit)] = 1;
                }
            }}));
    }

    for (auto &&task : tasks)
    {
        task.get();
    }
    uint64_t overall_count{0};
    for (uint64_t idx = 0; idx < bv_size; idx+=64)
    {
        overall_count += sdsl::bits::cnt(overall_content.get_int(idx));
    }
    std::cerr << "Overall" << '\t' << overall_count << std::endl;
}

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // check if file already exists or can be created
    if (!check_output_file(options.output_file))
        return 1;

    try
    {
        count_kmers(options);
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
