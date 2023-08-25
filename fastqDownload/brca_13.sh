for i in $(cat ~/tobedown.txt);
do
  ascp -QT -l 300m -P 33001 -i /home/zhepan/.aspera/connect/etc/aspera01.openssh aspera01@download.cncb.ac.cn:gsa-human/HRA002051/$i /home/zhepan/Project/PanCancerAtlas/data/BRCA_13/HRA002051/$i
done


#ascp -QT -l 300m -P 33001 -i /home/zhepan/.aspera/connect/etc/aspera01.openssh aspera01@download.cncb.ac.cn:gsa-human/HRA002051 /home/zhepan/Project/PanCancerAtlas/data/BRCA_13

#iptables -I INPUT -p udp --dport 33001 -j ACCEPT

#iptables -I OUTPUT -p udp --dport 33001 -j ACCEPT

