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

#include <atomic>
#include <set>

#include <seqan/arg_parse.h>
#include <seqan/binning_directory.h>

#include "helper.h"
#include "safequeue.hpp"

using namespace seqan;

struct ReadBatch
{
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
};

struct ReadResult
{
    CharString id;
    std::set<std::string> bins;
};

struct Options
{
    // CharString  contigs_dir;
    CharString  query_file;
    CharString  filter_file;
    CharString  output_file;

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

    addOption(parser, ArgParseOption("o", "output-file", "Specify an output filename for the results. \
                                     Default: search_results.txt", ArgParseOption::OUTPUT_FILE));

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
    // Parse contigs index prefix.
    getOptionValue(options.output_file, parser, "output-file");
    if (!isSet(parser, "output-file"))
    {
        options.output_file = CharString("search_results.txt");
    }

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
    std::array<uint32_t, 255> bin2file{
        0,0,0,0,
        1,1,1,1,1,1,
        2,2,2,2,2,2,
        3,3,3,3,3,
        4,4,4,4,
        5,5,5,5,5,5,
        6,6,6,6,6,
        7,7,7,7,7,
        8,8,8,8,8,
        9,9,9,9,9,9,
        10,10,10,10,10,
        11,11,11,11,11,
        12,12,12,12,12,12,12,
        13,13,13,13,13,13,13,13,13,13,
        14,14,14,14,14,14,14,14,14,14,
        15,15,15,15,15,15,15,15,15,15,15,
        16,16,16,16,16,16,16,16,16,16,16,16,
        17,17,17,17,17,17,17,17,17,17,17,
        18,18,18,18,
        19,19,19,19,
        20,20,20,20,
        21,21,21,
        22,22,22,
        23,23,
        24,24,24,24,
        25,25,25,25,25,25,
        26,26,26,26,
        27,27,27,
        28,28,28,28,
        29,29,29,29,
        30,30,30,30,30,30,30,
        31,31,31,31,31,31,
        32,32,32,32,32,32,
        33,33,33,33,33,
        34,34,34,34,
        35,35,35,
        36,36,36,
        37,37,37,
        38,38,38,38,
        39,39,39,
        40,40,40,
        41,41,41,41,
        42,42,42,42,
        43,43,
        44,44,44,44,
        45,45,45,45,45,45,45,45,
        46,46,46,46,
        47,
        48,48,48,48,
        49,49,49,49,49,49,49,49,49};
    std::array<std::string, 50> file2srr{
        "SRR1523653",
        "SRR1523654",
        "SRR1523655",
        "SRR1523656",
        "SRR1523657",
        "SRR1523658",
        "SRR1523659",
        "SRR1523661",
        "SRR1523662",
        "SRR1523663",
        "SRR1523664",
        "SRR1523665",
        "SRR1523666",
        "SRR2038259",
        "SRR2038310",
        "SRR2038322",
        "SRR2038440",
        "SRR2038441",
        "SRR5444611",
        "SRR5444613",
        "SRR5444615",
        "SRR5444617",
        "SRR5444619",
        "SRR5444621",
        "SRR5444623",
        "SRR5444625",
        "SRR5444643",
        "SRR5444645",
        "SRR5444647",
        "SRR5444649",
        "SRR5444651",
        "SRR5444653",
        "SRR5444655",
        "SRR5444657",
        "SRR5444661",
        "SRR5444665",
        "SRR5444669",
        "SRR5756304",
        "SRR5756312",
        "SRR5756317",
        "SRR5756320",
        "SRR5756324",
        "SRR5762372",
        "SRR5762373",
        "SRR5762374",
        "SRR5762375",
        "SRR5762376",
        "SRR5762377",
        "SRR5762378",
        "SRR5762379"
    };

    std::ofstream out(toCString(options.output_file));

    std::vector<std::future<void>> tasks;
    SafeQueue<ReadBatch> rbq;
    SafeQueue<ReadResult> rsq;

    bool finished_reading = false;

    tasks.emplace_back(std::async([=, &rbq, &options, &finished_reading] () {
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
            while (rbq.size() > options.threads * 3)
            {
                ;
            }
            StringSet<Dna5String> seq;
            StringSet<CharString> id;
            readRecords(id, seq, seq_file_in, 50000);
            rbq.push(ReadBatch{id, seq});
        }
        close(seq_file_in);
        finished_reading = true;
    }));
    std::atomic<uint32_t> finished_search{0};
    for (uint32_t task_number = 0; task_number < (options.threads > 2 ? options.threads - 2 : 1); ++task_number)
    {
        tasks.emplace_back(std::async([=, &rbq, &rsq, &filter, &file2srr, &bin2file, &finished_reading, &finished_search] () {
            while (true)
            {
                ReadBatch rb = rbq.pop();
                for (size_t i = 0; i < length(rb.ids); ++i)
                {
                    if(length(rb.seqs[i]) < getKmerSize(filter))
                        continue;
                    std::vector<bool> result = select(filter, rb.seqs[i], options.errors, options.penalty);
                    std::set<std::string> bins;
                    for (size_t j = 0; j < result.size(); ++j)
                    {
                        if (result[j])
                        {
                            bins.insert(file2srr[bin2file[j]]);
                        }
                    }
                    rsq.push(ReadResult{rb.ids[i], bins});
                }
                if (finished_reading && rbq.empty())
                {
                    finished_search++;
                    break;
                }
            }
        }));
    }

    tasks.emplace_back(std::async([=, &out, &rsq, &options, &finished_search] () {
        while (true)
        {
            ReadResult rs = rsq.pop();
            if (rs.id == "")
                continue;
            out << rs.id << '\n';
            if (!rs.bins.empty())
            {
                const auto separator = ",";
                const auto* sep = "";
                for(auto const & item : rs.bins) {
                    out << sep << item;
                    sep = separator;
                }
            }
            else
                out << "NA";
            out << std::endl;
            if (rsq.empty() && finished_search.load() == (options.threads > 2 ? options.threads - 2 : 1))
                break;
        }
    }));

    for (auto && task : tasks)
    {
         task.get();
    }
    out.close();
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
