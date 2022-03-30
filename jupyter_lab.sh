if  [ ! "$1" ] ;then
    thread_num=1
else
    thread_num=$1
fi

workdir=$(dirname $(readlink -f $0))

docker run -it --rm --user=$UID --gpus all -e OMP_NUM_THREADS=$thread_num -e PYTHONHASHSEED=0 -v $(pwd):/home/mydocker/$workdir -p 8888:8888 lyc1995/bioinfo:1.0.1_jupyter
