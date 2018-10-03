#!/bin/bash

# Install spark

INSTALL_DIR=${INSTALL_DIR:-/usr}
USER=`whoami`

SPARK=spark-${SPARK_VER:-2.1.1}-bin-hadoop${HADOOP_VER:-2.7}
SPARK_DIR=${INSTALL_DIR}/$SPARK
SPARK_LOCAL_DIR="/usr/local/spark"

download_spark() {
  wget -q http://d3kbcqa49mib13.cloudfront.net/$SPARK.tgz
  sudo tar -zxf $SPARK.tgz --directory $INSTALL_DIR &&
  sudo chown -R $USER:$USER $SPARK_DIR &&
  sudo ln -s $INSTALL_DIR/$SPARK $SPARK_LOCAL_DIR &&
  echo "download_spark successful"
}

configure_spark() {
  export SPARK_HOME=${SPARK_LOCAL_DIR} &&
  export PATH=${SPARK_LOCAL_DIR}/bin:${SPARK_LOCAL_DIR}/sbin:$PATH &&
  IP=`hostname -i` && MASTER=`hostname -s` &&
  sudo echo "SPARK_MASTER_HOST=$IP" > ${SPARK_LOCAL_DIR}/conf/spark-env.sh &&
  sudo echo "SPARK_LOCAL_IP=$IP" >> ${SPARK_LOCAL_DIR}/conf/spark-env.sh &&
  sudo echo "localhost" > ${SPARK_LOCAL_DIR}/conf/slaves &&
  echo "$IP $MASTER" | sudo tee --append /etc/hosts &&
  echo "configure_spark successful"
}

install_spark() {
  download_spark &&
  configure_spark &&
  echo "Install Spark successful"
}

install_spark
