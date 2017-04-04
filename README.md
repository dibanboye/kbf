### Project Proposal 
Alan Kuhnle, Victoria Crawford, Yuanpu Xie, Zizhao Zhang, Ke Bo

Compact data structures for the representation of de Bruijn graphs are important for genome sequence assembly. A data structure supporting efficient and exact membership queries for fully dynamic de Bruijn graphs was recently proposed in Belazzougui et al. (2016). For our project, we will implement this data structure and evaluate its practical performance as compared to popular implementations based upon bloom filters when used for de Bruijn construction from sequence data.


##Instruction to run:
- Go to ./cpp-src
- Compile
    ```
    sh compile.sh main
    ```
- Test
    ```
    ./main data/yeast.fasta 10 | tee log/terminal_log.txt
    ```

log.txt is written in log/.

terminal_log.txt save all system output for python parsing.
