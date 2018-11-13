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

#include <seqan/arg_parse.h>
#include <seqan/binning_directory.h>

#include "helper.h"

using namespace seqan;

struct Options
{
    // CharString  contigs_dir;
    CharString  query_file;
    CharString  filter_file;

    uint32_t    errors;
    uint32_t    penalty;
    // uint32_t    kmer_size;
    uint32_t    window_size;
    // uint32_t    number_of_bins;
    // uint64_t    size_of_ibf;
    // uint32_t    number_of_hashes;
    unsigned    threads;

    Options():
        errors(0),
        penalty(0),
        // kmer_size(19),
        window_size(24),
        // number_of_bins(64),
        // size_of_ibf(16_g),
        // number_of_hashes(3),
        threads(1) {}
};

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "SRA_search build prototype");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "QUERY FILE"));
    setHelpText(parser, 0, "A file containing the reads to query");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "IBF FILE"));
    setHelpText(parser, 0, "A file containing the IBF to query");


    addSection(parser, "Query Options");

    // addOption(parser, ArgParseOption("o", "output-file", "Specify an output filename for the filter. \
    //                                  Default: use the directory name of reference genomes.", ArgParseOption::OUTPUT_FILE));
    // setValidValues(parser, "output-file", "filter");

    // addOption(parser, ArgParseOption("b", "number-of-bins", "The number of bins",
    //                                  ArgParseOption::INTEGER));

    // setMinValue(parser, "number-of-bins", "1");
    // setMaxValue(parser, "number-of-bins", "4194300");

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threads);

    addOption(parser, ArgParseOption("e", "errors", "Maximum number of errors to allow.", ArgParseOption::INTEGER));
    setMinValue(parser, "errors", "0");
    setMaxValue(parser, "errors", "10");
    setDefaultValue(parser, "errors", options.errors);

    addOption(parser, ArgParseOption("p", "penalty", "Correctional value for threshold calculation", ArgParseOption::INTEGER));
    setMinValue(parser, "penalty", "0");
    setMaxValue(parser, "penalty", "10");
    setDefaultValue(parser, "penalty", options.penalty);

    // addOption(parser, ArgParseOption("k", "kmer-size", "The size of kmers for the IBF",
    //                                  ArgParseOption::INTEGER));
    // setMinValue(parser, "kmer-size", "14");
    // setMaxValue(parser, "kmer-size", "32");

    addOption(parser, ArgParseOption("w", "window-size", "The size of the window for the IBF",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "window-size", "14");

    // addOption(parser, ArgParseOption("nh", "num-hash", "Specify the number of hash functions to use for the bloom filter.", ArgParseOption::INTEGER));
    // setMinValue(parser, "num-hash", "2");
    // setMaxValue(parser, "num-hash", "5");
    // setDefaultValue(parser, "num-hash", options.number_of_hashes);

    // addOption(parser, ArgParseOption("bs", "bloom-size",
    //         "The size of bloom filter suffixed by either M or G for megabytes or gigabytes respectively.",
    //         ArgParseOption::STRING));
    // setDefaultValue(parser, "bloom-size", "1G");
}

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.query_file, parser, 0);
    getArgumentValue(options.filter_file, parser, 1);

    // // Parse contigs input file.
    // getArgumentValue(options.contigs_dir, parser, 0);
    //
    // // Append trailing slash if it doesn't exist.
    // append_trailing_slash(options.contigs_dir);
    //
    // // Parse contigs index prefix.
    // getOptionValue(options.filter_file, parser, "output-file");
    // if (!isSet(parser, "output-file"))
    // {
    //     options.filter_file = trimExtension(options.contigs_dir);
    //     append(options.filter_file, "bloom.filter");
    // }

    if (isSet(parser, "errors")) getOptionValue(options.errors, parser, "errors");
    if (isSet(parser, "penalty")) getOptionValue(options.penalty, parser, "penalty");
    // if (isSet(parser, "number-of-bins")) getOptionValue(options.number_of_bins, parser, "number-of-bins");
    // if (isSet(parser, "kmer-size")) getOptionValue(options.kmer_size, parser, "kmer-size");
    if (isSet(parser, "window-size")) getOptionValue(options.window_size, parser, "window-size");
    if (isSet(parser, "threads")) getOptionValue(options.threads, parser, "threads");
    // if (isSet(parser, "num-hash")) getOptionValue(options.number_of_hashes, parser, "num-hash");

    // std::string ibf_size;
    // if (getOptionValue(ibf_size, parser, "bloom-size"))
    // {
    //     uint64_t base = std::stoi(ibf_size);
    //     switch (ibf_size.at(ibf_size.size()-1))
    //     {
    //         case 'G': case 'g':
    //             options.size_of_ibf =  base * 8*1024*1024*1024;
    //             break;
    //         case 'M': case 'm':
    //             options.size_of_ibf =  base * 8*1024*1024;
    //             break;
    //         default:
    //             std::cerr <<"[ERROR] invalid --bloom-size (-bs) parameter provided. (eg 256M, 1g)" << std::endl;
    //             exit(1);
    //     }
    //
    // }
    return ArgumentParser::PARSE_OK;
}

template <typename TFilter>
inline void search_filter(Options & options, TFilter & filter)
{
    Dna5String seq;
    CharString id;
    SeqFileIn seq_file_in;
    if (!open(seq_file_in, toCString(options.query_file)))
    {
        CharString msg = "Unable to open contigs file: ";
        append(msg, CharString(options.query_file));
        std::cerr << msg << std::endl;
        throw toCString(msg);
    }
    while(!atEnd(seq_file_in))
    {
        readRecord(id, seq, seq_file_in);
        if(length(seq) < getKmerSize(filter))
            continue;
        std::vector<uint64_t> result = count(filter, seq);
        std::cerr << "Read " << id << " counts in bins:\n";
        for (size_t i = 0; i < result.size(); ++i)
        {
            std::cerr << i << ' ' << result[i] << " || ";
        }
        std::cerr << std::endl;
        // bool found{false};
        // std::vector<bool> result = select(filter, seq, options.errors, options.penalty);
        // std::cerr << "Read " << id << " found in bins:\n";
        // for (size_t i = 0; i < result.size(); ++i)
        // {
        //     if (result[i])
        //     {
        //         found = true;
        //         std::cerr << i << ' ';
        //     }
        // }
        // if (!found)
        //     std::cerr << "None";
        // std::cerr << std::endl;
    }
    // std::string com_ext = common_ext(options.contigs_dir, options.number_of_bins);
    //
    // uint32_t batch_size = options.number_of_bins/options.threads;
    // if(batch_size * options.threads < options.number_of_bins) ++batch_size;
    //
    // std::vector<std::future<void>> tasks;
    //
    // for (uint32_t task_number = 0; task_number < options.threads; ++task_number)
    // {
    //     tasks.emplace_back(std::async([=, &filter] {
    //         for (uint32_t bin_number = task_number*batch_size;
    //             bin_number < options.number_of_bins && bin_number < (task_number +1) * batch_size;
    //             ++bin_number)
    //         {
    //             CharString seq_file_path;
    //             append_file_name(seq_file_path, options.contigs_dir, bin_number);
    //             append(seq_file_path, com_ext);
    //
    //             // read everything as CharString to avoid impure sequences crashing the program
    //             CharString id, seq;
    //             SeqFileIn seq_file_in;
    //             if (!open(seq_file_in, toCString(seq_file_path)))
    //             {
    //                 CharString msg = "Unable to open contigs file: ";
    //                 append(msg, CharString(seq_file_path));
    //                 std::cerr << msg << std::endl;
    //                 throw toCString(msg);
    //             }
    //             while(!atEnd(seq_file_in))
    //             {
    //                 readRecord(id, seq, seq_file_in);
    //                 if(length(seq) < options.kmer_size)
    //                     continue;
    //                 insertKmer(filter, seq, bin_number);
    //             }
    //         }}));
    // }
    //
    // for (auto &&task : tasks)
    // {
    //     task.get();
    // }
    // store(filter, toCString(options.filter_file));
}

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // // check if file already exists or can be created
    // if (!check_output_file(options.filter_file))
    //     return 1;

    try
    {
        typedef BDConfig<Dna5, Minimizer<19, 24>, Uncompressed> Config;
        BinningDirectory<InterleavedBloomFilter, Config> const filter(options.filter_file, options.window_size);
        search_filter(options, filter);
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
