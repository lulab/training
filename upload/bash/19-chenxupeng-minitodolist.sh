#!/bin/bash
# a bash script to create a mini todo list
addTask()
{
if [ ! -f ~/.todo/list.txt ]; then
if [ ! -d ~/.todo ]; then mkdir ~/.todo; fi
touch ~/.todo/list.txt
fi
echo "$1       ;      $(date)" >> ~/.todo/list.txt
}

checkInternet()
{
httpGet github.com > /dev/null 2>&1 || { echo "Error: no active internet connection" >&2; return 1; } # query github with a get request
}

removeTask()
{
## Check for valid task numbers (valid characters)
if [ -f ~/.todo/temp.txt ];then rm -f ~/.todo/temp.txt;fi
touch ~/.todo/temp.txt
for taskToRemove in "$@";do
oldTaskNumber=$taskToRemove
taskNumber=$( echo $taskToRemove | grep -Eo "[0-9]*" )
if [[ $taskNumber == "" || $oldTaskNumber != $taskNumber ]]; then echo "Error: $oldTaskNumber is not a valid task number!" && return 1; fi
done
count="0"
IFS=$'\n'       # make newlines the only separator

## Removing the task (only don't add to temp if we should remove it)
for task in $(cat ~/.todo/list.txt); do
removeIt="false"
for taskToRemove in "$@";do
if [[ $(($count + 1)) == "$taskToRemove" ]]; then
removeIt="true"
break
fi
done
if ! $removeIt ;then echo "$task" >> ~/.todo/temp.txt;fi
count=$(( $count + 1 ))
done
rm -f ~/.todo/list.txt
cp  ~/.todo/temp.txt ~/.todo/list.txt
rm -f ~/.todo/temp.txt

##Checking if the task exists
for taskToRemove in "$@" ;do
if [ $count -lt $taskToRemove ]; then
echo "Error: task number $taskToRemove does not exist!"
else
echo "Sucessfully removed task number $taskToRemove"
fi
done
}

getTasks()
{
if [ -f ~/.todo/list.txt ]; then
checkEmpty=$(cat ~/.todo/list.txt)
if [[ $checkEmpty == "" ]]; then
echo "No tasks found"
else
count="1"
IFS=$'\n'       # make newlines the only separator
for task in $(cat ~/.todo/list.txt); do
tempTask=$count
if [ $count -lt 10 ]; then tempTask="0$count"; fi
echo "$tempTask). $task"  >> ~/.todo/getTemp.txt
count=$(( $count + 1 ))
done
cat ~/.todo/getTemp.txt | column -t -s ";"
rm -f ~/.todo/getTemp.txt
fi
else
echo "No tasks found"
fi
}

usage()
{
cat <<EOF
-a  Add the following task
-g  Get the current tasks
-r  Remove by task numbers
EOF
}

while getopts "cr:a:guvh" opt; do
case "$opt" in
g)  if [[ $flag == "" ]]; then
flag="get"
else
echo "Error: all flags are mutually exclusive"
exit 1
fi
;;
r)  if [[ $flag == "" ]]; then
flag="remove"
else
echo "Error: all flags are mutually exclusive"
exit 1
fi
;;
a)  if [[ $flag == "" ]]; then
flag="add"
else
echo "Error: all flags are mutually exclusive"
exit 1
fi
;;
esac
done

if [[ $# == "0" ]]; then
usage
elif [[ $# == "1" ]]; then
if [[ $flag == "get" || $1 == "list" || $1 == "get" ]]; then getTasks || exit 1
else { echo "Error: the argument $1 is not valid"; exit 1; }; fi
else
if [[ $flag == "add" || $1 == "add" ]]; then addTask "${*:2}" && getTasks || exit 1
elif [[ $flag == "remove" || $1 == "remove" ]]; then removeTask ${*:2} && getTasks || exit 1
else { usage; exit 1; }; fi
fi
