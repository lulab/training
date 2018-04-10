#!/bin/bash

linenum=`df | egrep "\<[[:digit:]]{1,3}%.*" | wc -l`
test1=`df | egrep -o "\<[[:digit:]]{1,3}%.*" | tr '%' ' ' | tr -s ' ' |sort -nr | head -n1 | cut -d' ' -f1`
test2=`df -i | egrep -o "\<[[:digit:]]{1,3}%.*" | tr '%' ' ' | tr -s ' ' |sort -nr | head -n1 | cut -d' ' -f1`

num=1

if [ "$test1" -ge "80" ];then
    while [ "$num" -le "$linenum" ]
    do
        testnum1=`df | egrep -o "\<[[:digit:]]{1,3}%.*" | tr '%' ' ' | tr -s ' ' |sort -nr | head -n"$num" | tail -n1 | head -n1 | cut -d' ' -f1`
        #提取出分区使用量的数值用于后面的比较
        testname1=`df | egrep -o "\<[[:digit:]]{1,3}%.*" | tr '%' ' ' | tr -s ' ' |sort -nr | head -n"$num" | tail -n1 | head -n1 | cut -d' ' -f2`
        #提取分区路径
        [ "$testnum1" -ge "80" ] && wall "$testname1"' will be full'
        #比较数值的和输出的部分
        num=$[num+1]
    done

num=1

elif [ "$test2" -ge "80" ];then
    while [ "$num" -le "$linenum" ]
    do
        testnum2=`df -i | egrep -o "\<[[:digit:]]{1,3}%.*" | tr '%' ' ' | tr -s ' ' |sort -nr | head -n"$num" | tail -n1 | head -n1 | cut -d' ' -f1`
        testname2=`df -i | egrep -o "\<[[:digit:]]{1,3}%.*" | tr '%' ' ' | tr -s ' ' |sort -nr | head -n"$num" | tail -n1 | head -n1 | cut -d' ' -f2`
        [ "$testnum2" -ge "80" ] && wall "$testname2"' Inode will be full'
        num=$[num+1]
    done

else
    echo "All filesystems are safe"
fi

unset linenum
unset test1
unset test2
unset num
unset testnum1
unset testname1
unset testnum2
unset testname2
exit
