/**
 * The MIT License (MIT)
 * Copyright (c) 2019 Omics Data Automation, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of 
 * this software and associated documentation files (the "Software"), to deal in 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
 * the Software, and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


// Amortize the cost of building catch2.hpp for
// every unit test by linking in ctest_main.o

//See https://github.com/catchorg/Catch2/blob/master/docs/own-main.md#adding-your-own-command-line-options
//to pass run time arguments to tests

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

std::string g_query_json_file = "";
std::string g_loader_json_file = "";
std::string g_golden_output_file = "";
bool g_skip_GT_matching = false;

int main( int argc, char* argv[] )
{
  Catch::Session session; // There must be exactly one instance

  // Build a new parser on top of Catch's
  using namespace Catch::clara;
  auto cli 
    = session.cli() // Get Catch's composite command line parser
    | Opt( g_query_json_file, "Query JSON file" ) // bind variable to a new option, with a hint string
     ["--query-json-file"]    // the option names it will respond to
    ("Query JSON file")        // description string for the help output
    | Opt(g_loader_json_file, "Loader JSON file")
     ["--loader-json-file"]
     ("Loader JSON file")
    | Opt(g_golden_output_file, "Golden output file")
     ["--golden-output-file"]
     ("Golden output file")
    | Opt(g_skip_GT_matching, "Skip GT matching")
     ["--skip-GT-matching"]
     ("Skip GT matching")
    ;

  // Now pass the new composite back to Catch so it uses that
  session.cli( cli ); 

  // Let Catch (using Clara) parse the command line
  int returnCode = session.applyCommandLine( argc, argv );
  if( returnCode != 0 ) // Indicates a command line error
    return returnCode;

  return session.run();
}
