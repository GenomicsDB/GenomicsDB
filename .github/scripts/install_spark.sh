#!/bin/bash

# Install spark

INSTALL_DIR=${INSTALL_DIR:-/usr}
USER=`whoami`

SPARK_VER=${SPARK_VER:-2.4.0}
SPARK=spark-$SPARK_VER-bin-hadoop${SPARK_HADOOP_VER:-2.7}
SPARK_DIR=${INSTALL_DIR}/$SPARK
SPARK_LOCAL_DIR="/usr/local/spark"
SPARK_ENV=${SPARK_ENV:-$HOME/spark_env.sh}

# retry logic from: https://docs.microsoft.com/en-us/azure/hdinsight/hdinsight-hadoop-script-actions-linux
MAXATTEMPTS=3
retry() {
    local -r CMD="$@"
    local -i ATTEMPTNUM=1
    local -i RETRYINTERVAL=2

    until $CMD
    do
        if (( ATTEMPTNUM == MAXATTEMPTS ))
        then
                echo "Attempt $ATTEMPTNUM failed. no more attempts left."
                return 1
        else
                echo "Attempt $ATTEMPTNUM failed! Retrying in $RETRYINTERVAL seconds..."
                sleep $(( RETRYINTERVAL ))
                ATTEMPTNUM=$ATTEMPTNUM+1
        fi
    done
}

download_spark() {
  retry wget -nv --trust-server-names "https://archive.apache.org/dist/spark/spark-$SPARK_VER/$SPARK.tgz"
  sudo tar -zxf $SPARK.tgz --directory $INSTALL_DIR &&
  sudo chown -R $USER:$USER $SPARK_DIR &&
  sudo ln -s $INSTALL_DIR/$SPARK $SPARK_LOCAL_DIR &&
  echo "download_spark successful"
}

setup_spark_env() {
  echo "export SPARK_HOME=${SPARK_LOCAL_DIR}" >> $SPARK_ENV &&
  echo "export PATH=${SPARK_LOCAL_DIR}/bin:${SPARK_LOCAL_DIR}/sbin:$PATH" >> $SPARK_ENV &&
  echo "export CLASSPATH=${SPARK_LOCAL_DIR}/jars/gcs-connector-latest-hadoop2.jar:$CLASSPATH" >> $SPARK_ENV &&
  source $SPARK_ENV
}

configure_spark() {
  echo "export SPARK_HOME=${SPARK_LOCAL_DIR}" >> $SPARK_ENV &&
  echo "export PATH=${SPARK_LOCAL_DIR}/bin:${SPARK_LOCAL_DIR}/sbin:$PATH" >> $SPARK_ENV &&
  echo "export CLASSPATH=${SPARK_LOCAL_DIR}/jars/gcs-connector-latest-hadoop2.jar:$CLASSPATH" >> $SPARK_ENV &&
  source $SPARK_ENV &&
  IP=`hostname -i` && MASTER=`hostname -s` &&
  sudo echo "SPARK_MASTER_HOST=$IP" > ${SPARK_LOCAL_DIR}/conf/spark-env.sh &&
  sudo echo "SPARK_LOCAL_IP=$IP" >> ${SPARK_LOCAL_DIR}/conf/spark-env.sh &&
  sudo echo "localhost" > ${SPARK_LOCAL_DIR}/conf/slaves &&
  sudo cp ${SPARK_LOCAL_DIR}/conf/log4j.properties.template ${SPARK_LOCAL_DIR}/conf/log4j.properties &&
  echo "$IP $MASTER" | sudo tee --append /etc/hosts &&
  echo "configure_spark successful"
}

install_spark() {
  if [[ ! -f ${SPARK_DIR}/conf/slaves ]]; then
    echo "Installing Spark..."
    download_spark &&
    configure_spark &&
    echo "Install Spark successful"
  else
    echo "Found cached Spark install" 
  fi
}

install_spark &&
if [[ ! -L ${SPARK_LOCAL_DIR} ]]; then sudo ln -s $INSTALL_DIR/$SPARK $SPARK_LOCAL_DIR; fi &&
setup_spark_env &&
${SPARK_LOCAL_DIR}/sbin/start-master.sh &&
echo "Started spark"
