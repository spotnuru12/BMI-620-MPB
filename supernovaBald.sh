mkdir -p ../supergnova_out

declare -A NSIZE=(
  ["Insomnia"]=386533
  ["Cannabis"]=357806
  ["Income"]=286301
  ["Chronotype"]=449734
  ["ADHD"]=225534
)
# output 5 seperate text files -> stitch together at the end
for TRAIT in Insomnia Cannabis Income Chronotype ADHD; do
  echo ""
  echo "=== Baldness vs $TRAIT ==="
  python3 supergnova.py \
    "../Baldness_sumstats.gz" \
    "../${TRAIT}_sumstats.gz" \
    --N1 205327 \
    --N2 ${NSIZE[$TRAIT]} \
    --bfile data/bfiles/eur_chr@_SNPmaf5 \
    --partition data/partition/eur_chr@.bed \
    --thread 4 \
    --out ../supergnova_out/Baldness_vs_${TRAIT}.txt
done

echo ""
echo "=== All done! Files: ==="
ls -la ../supergnova_out/
