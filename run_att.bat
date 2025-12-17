@echo off
REM Windows batch script to run process_att.m

set "THRESHOLDS=0.988 0.993 0.996"
set "WAIT_FOR_COMPLETION=0"

if "%WAIT_FOR_COMPLETION%"=="1" (
    set "CMD_MAT=matlab -nodisplay -nosplash -wait -r"
) else (
    set "CMD_MAT=matlab -nodisplay -nosplash -r"
)

for %%t in (%THRESHOLDS%) do (
    echo Running with threshold: %%t
    %CMD_MAT% "try, process_att(%%t); catch e, disp(e.message); end; exit;"
)
