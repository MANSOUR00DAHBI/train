# makeFile main
clear 
gcc -S  main.c 
gcc -c  main.s 
gcc -Wall -Wextra -o Setup main.o
rm *.s *.o 
./Setup 