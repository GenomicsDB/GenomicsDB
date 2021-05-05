#!/bin/bash

# Install spark

INSTALL_DIR=${INSTALL_DIR:-/usr}
USER=`whoami`

SPARK_VER=${SPARK_VER:-3.0.1}
SPARK=spark-$SPARK_VER-bin-hadoop${HADOOP_VER:-2.7}
SPARK_DIR=${INSTALL_DIR}/$SPARK
SPARK_LOCAL_DIR="/usr/local/spark"

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
  if [[ ! -f $CACHE_DIR/$SPARK.tgz ]]; then
    retry wget -nv --trust-server-names "https://archive.apache.org/dist/spark/spark-$SPARK_VER/$SPARK.tgz" -O $CACHE_DIR/$SPARK.tgz
  fi
  sudo tar -zxf $CACHE_DIR/$SPARK.tgz --directory $INSTALL_DIR &&
  sudo chown -R $USER:$USER $SPARK_DIR &&
  sudo ln -s $INSTALL_DIR/$SPARK $SPARK_LOCAL_DIR &&
  get_gcs_connector &&
  echo "download_spark successful"
}

get_gcs_connector() {
  if [[ $INSTALL_TYPE == gcs ]]; then
    mv "$TRAVIS_BUILD_DIR/gcs-connector-latest-hadoop2.jar" "$SPARK_LOCAL_DIR/jars"
  fi
}

configure_spark() {
  export SPARK_HOME=${SPARK_LOCAL_DIR} &&
  export PATH=${SPARK_LOCAL_DIR}/bin:${SPARK_LOCAL_DIR}/sbin:$PATH &&
  export CLASSPATH=${SPARK_LOCAL_DIR}/jars/gcs-connector-latest-hadoop2.jar:$CLASSPATH &&
  IP=`hostname -i` && MASTER=`hostname -s` &&
  sudo echo "SPARK_MASTER_HOST=$IP" > ${SPARK_LOCAL_DIR}/conf/spark-env.sh &&
  sudo echo "SPARK_LOCAL_IP=$IP" >> ${SPARK_LOCAL_DIR}/conf/spark-env.sh &&
  sudo echo "localhost" > ${SPARK_LOCAL_DIR}/conf/slaves &&
  sudo cp ${SPARK_LOCAL_DIR}/conf/log4j.properties.template ${SPARK_LOCAL_DIR}/conf/log4j.properties &&
  echo "$IP $MASTER" | sudo tee --append /etc/hosts &&
  echo "configure_spark successful"
}

install_spark() {
  echo "Installing Spark..."
  download_spark &&
  configure_spark &&
  echo "Install Spark successful"
}

install_spark
