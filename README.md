**1.Run FEMU Code**


**安装依赖编译**

```
cd femu
mkdir build-femu
cd build-femu
cp ../femu-scripts/femu-copy-scripts.sh .
./femu-copy-scripts.sh .
# only Debian/Ubuntu based distributions supported
sudo ./pkgdep.sh

./femu-compile.sh
```

**vim femu-compile.sh run-blackbox.sh** 

```
# femu-compile.sh

#!/bin/bash
  
NRCPUS="$(cat /proc/cpuinfo | grep "vendor_id" | wc -l)"

make clean
# --disable-werror --extra-cflags=-w --disable-git-update
../configure --disable-werror --enable-kvm --target-list=x86_64-softmmu
make -j $NRCPUS

echo ""
echo "===> FEMU compilation done ..."
echo ""
exit


# run-blackbox.sh
if [[ ! -e "$OSIMGF" ]]; then
        echo ""
        echo "VM disk image couldn't be found ..."
        echo "Please prepare a usable VM image and place it as $OSIMGF"
        echo "Once VM disk image is ready, please rerun this script again"
        echo ""
        exit
fi

sudo x86_64-softmmu/qemu-system-x86_64 \
    -name "FEMU-BBSSD-VM" \
    -enable-kvm \
    -cpu host \
    -smp 4 \
    -m 12G \
    -device virtio-scsi-pci,id=scsi0 \
    -device scsi-hd,drive=hd0 \
    -drive file=$OSIMGF,if=none,aio=native,cache=none,format=qcow2,id=hd0 \
    -device femu,devsz_mb=16384,femu_mode=1 \
    -net user,hostfwd=tcp::8083-:22 \
    -net nic,model=virtio \
    -nographic \
    -qmp unix:./qmp-sock,server,nowait 2>&1 | tee log
```





**启动**

```
./run-blackbox.sh

# if you want to connect in another terminal
ssh -p 8083 femu@localhost
# femu 的密码：femu
```



**mkfs & mount**

```
mkdir test
sudo mkfs.ext4 /dev/nvme0n1
sudo mount /dev/nvme0n1 test/

# use lsblk & df -h to check whether is mouted successfully


```



**write fio scripts**

```
[global]
ioengine=sync
direct=1
rw=randwrite
bs=4k
size=21M
numjobs=1

[job1]
filename=/dev/nvme0n1
do_verify=0
write_iolog=/home/femu/trace/FIU/homes/homes-110108-112108.1.blkparse
```





**nvme-cli  && log**

```
# you can use nvme-cli to change ssd_write interface
sudo nvme admin-passthru /dev/nvme0n1 --opcode=0xef --cdw10=8


# Log file path
# In your host not in femu
cd ~/Femu-caftl/build-femu/
cat log
```


