THRESHOLDS="0.988 0.993 0.996"
WAIT_FOR_COMPLETION=0

if [ "$WAIT_FOR_COMPLETION" -eq "1" ]; then
    cmd_mat_ng="matlab -nodisplay -nosplash -wait -r"
else
    cmd_mat_ng="matlab -nodisplay -nosplash -r"
fi

for t in $THRESHOLDS; do
    echo "Running with threshold: $t"
    $cmd_mat_ng  "BIN_DIR = '/home1/cyi/workplace/wks_on_nobkup/povray/isosurface/BIN'; addpath(genpath(BIN_DIR));process_att($t); exit;"
done

