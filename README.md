# CPU Cache Simulation
Performs a simulation of a specified algorithm on an emulated cache and prints statistics.

## Example Usage
```
./a.out -c 262144 -b 32
INPUTS=====================================
Ram Size =                   5529600 bytes
Cache Size =                 262144 bytes
Block Size =                 32 bytes
Total Blocks in Cache =      8192
Associativity =              2
Number of Sets =             4096
Replacement Policy =         LRU
Algorithm =                  blocked mxm
MXM Blocking Factor =        32
Matrix or Vector Dimension = 480
RESULTS=====================================
Instruction Count =          442598400
Read Hits =                  219589197
Read Misses =                1594803
Read Miss Rate =             0.721030 %
```

## Parameters
| Parameter                 | Description       | Default   |	
| :------------------------ |:-------------:| :-------------|
| -c | cache size (in bytes) | 65536
| -b | cache data block size (in bytes) | 64
| -n | cache associativity | 2
| -r | replacement policy (0: random, 1: FIFO, 2: LRU) | 2
| -a | algorithm (0: vector multiply, 1: matrix multiply, 2: blocked matrix multiply) | 2
| -d | matrix dimension | 480
| -f | blocking factor | 32
| -p | print matrix | 0

## Build

g++ *.cpp -std=c++11

## Analysis - Associativity

| Cache Associativity                 | Cache Size (Bytes) | Instructions       | Read Hits   |	Read Misses | Read Miss % |
| :------------------------ |:-------------:|:-------------:|:-------------:|:-------------:|:-------------|
|1|262144|442598400|217085427|4098573|1.85
|2|262144|442598400|220385397|798603|0.36
|8|262144|442598400|220336427|847573|0.38

## Analysis - Memory Block Size

| Memory Block Size (bytes)               | Cache Size (Bytes) | Instructions       | Read Hits   |	Read Misses | Read Miss % |
| :------------------------ |:-------------:|:-------------:|:-------------:|:-------------:|:-------------|
|8|262144|442598400|214812109|6371891|2.88
|32|262144|442598400|219589213|1594787|0.72
|1024|262144|442598400|109604540|111579460|50.44

## Analysis - Total Cache Size

| Cache Size (bytes)               | Instructions       | Read Hits   |	Read Misses | Read Miss % |
| :------------------------ |:-------------:|:-------------:|:-------------:|:-------------|
|4096|442598400|106488000|114696000|51.86
|8192|442598400|106522170|114661830|51.84
|16384|442598400|106550685|114633315|51.83
|32768|442598400|106614049|114569951|51.80
|65536|442598400|106746532|114437468|51.74
|131072|442598400|210635318|114437468|4.77
|262144|442598400|220385397|798603|0.36
|524288|442598400|220633778|550222|0.24

