#! /bin/bash
modes=("rows" "clms" "blocks")
for mode in ${modes[@]}
do
    mpicc main.c utils.c parallel-product-${mode}.c -o main
    for ((processes=1; processes<=8; processes*=2))
    do
        for ((data=512; data<=8192; data*=2))
        do
            mpiexec -n=${processes} ./main ${data} ${mode}
        done
    done
done
