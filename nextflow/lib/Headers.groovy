/*
 * This file holds several functions used to render the nf-core ANSI header.
 */


// MIT License
// Copyright (c) 2018 nf-core
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

class Headers {

    private static Map log_colours(Boolean monochrome_logs) {
        Map colorcodes = [:]
        colorcodes['reset']       = monochrome_logs ? '' : "\033[0m"
        colorcodes['dim']         = monochrome_logs ? '' : "\033[2m"
        colorcodes['black']       = monochrome_logs ? '' : "\033[0;30m"
        colorcodes['green']       = monochrome_logs ? '' : "\033[0;32m"
        colorcodes['yellow']      = monochrome_logs ? '' :  "\033[0;33m"
        colorcodes['yellow_bold'] = monochrome_logs ? '' : "\033[1;93m"
        colorcodes['blue']        = monochrome_logs ? '' : "\033[0;34m"
        colorcodes['purple']      = monochrome_logs ? '' : "\033[0;35m"
        colorcodes['cyan']        = monochrome_logs ? '' : "\033[0;36m"
        colorcodes['white']       = monochrome_logs ? '' : "\033[0;37m"
        colorcodes['red']         = monochrome_logs ? '' : "\033[1;91m"
        return colorcodes
    }

    static String dashed_line(monochrome_logs) {
        Map colors = log_colours(monochrome_logs)
        return "-${colors.dim}----------------------------------------------------${colors.reset}-"
    }

    static String nf_core(workflow, monochrome_logs) {
        Map colors = log_colours(monochrome_logs)
        String.format(
            """\n
            ${dashed_line(monochrome_logs)}
                                                    ${colors.green},--.${colors.black}/${colors.green},-.${colors.reset}
            ${colors.blue}        ___     __   __   __   ___     ${colors.green}/,-._.--~\'${colors.reset}
            ${colors.blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${colors.yellow}}  {${colors.reset}
            ${colors.blue}  | \\| |       \\__, \\__/ |  \\ |___     ${colors.green}\\`-._,-`-,${colors.reset}
                                                    ${colors.green}`._,._,\'${colors.reset}
            ${colors.purple}  ${workflow.manifest.name} v${workflow.manifest.version}${colors.reset}
            ${dashed_line(monochrome_logs)}
            """.stripIndent()
        )
    }
}
