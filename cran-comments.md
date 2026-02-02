## R CMD check results

0 errors | 0 warnings | 1 note

* This is a update release.

# Test results

There were no ERRORs or WARNINGs on R-hub Ubuntu, macos, and Windows, and local macOS.

There is one one NOTE which is explained below:

"A complete check needs the 'checkbashisms' script.""

Explanation: This warning appeared during local checks on macOS because the checkbashisms utility is not present on the local system. However, the package was tested on R-hub (Debian Linux), where the script is available, and it passed without any bashism issues. The configure script follows standard sh syntax.

