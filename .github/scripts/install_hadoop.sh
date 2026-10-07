#!/bin/bash

# Install hadoop
# Installation relies on finding JAVA_HOME@/usr/java/latest as a prerequisite

INSTALL_DIR=${INSTALL_DIR:-/usr}
USER=`whoami`

HADOOP=hadoop-${HADOOP_VER:-3.3.5}
HADOOP_DIR=${INSTALL_DIR}/$HADOOP
HADOOP_ENV=${HADOOP_ENV:-$HOME/hadoop_env.sh}
JAVA_VER=${JAVA_VER:-17}

install_prereqs() {
  # Java install in workflow yaml
  if [[ -f /usr/java/latest ]]; then
    echo "/usr/java/latest found"
    sudo rm /usr/java/latest
  fi
  if [[ ! -z $JAVA_HOME ]]; then
    sudo mkdir -p /usr/java
    sudo ln -s $JAVA_HOME /usr/java/latest
  else
    sudo ln -s /usr/lib/jvm/java-1.${JAVA_VER}.0-openjdk-amd64/ /usr/java/latest
  fi
  echo "install_prereqs successful"
}

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

# Download $1, a path under Apache's dist/, to $2. dlcdn.apache.org is fast but only carries current
# releases; archive.apache.org has every release but is slow and can stall mid-download, so --timeout
# turns a stall into a failed attempt that retry can repeat.
download_apache() {
  wget -nv --timeout=60 -O $2 https://dlcdn.apache.org/$1 ||
    retry wget -nv --timeout=60 -O $2 https://archive.apache.org/dist/$1
}

download_hadoop() {
  download_apache hadoop/common/$HADOOP/$HADOOP.tar.gz $HADOOP.tar.gz &&
  tar -xzf $HADOOP.tar.gz --directory $INSTALL_DIR &&
  echo "download_hadoop successful"
}

configure_passphraseless_ssh() {
  sudo apt update; sudo apt -y install openssh-server
  cat > sshd_config << EOF
          SyslogFacility AUTHPRIV
          PermitRootLogin yes
          AuthorizedKeysFile	.ssh/authorized_keys
          PasswordAuthentication yes
          ChallengeResponseAuthentication no
          UsePAM yes
          UseDNS no
          X11Forwarding no
          PrintMotd no
EOF
  sudo mv sshd_config /etc/ssh/sshd_config &&
  sudo systemctl restart ssh &&
  rm ~/.ssh/id_rsa 2> /dev/null
  ssh-keygen -q -t rsa -b 4096 -N '' -f ~/.ssh/id_rsa &&
  cat ~/.ssh/id_rsa.pub | tee -a ~/.ssh/authorized_keys &&
  chmod 600 ~/.ssh/authorized_keys &&
  chmod 700 ~/.ssh &&
  sudo chmod -c 0755 ~/ &&
  echo "configure_passphraseless_ssh successful"
}

configure_hadoop() {
  configure_passphraseless_ssh &&
  cp -fr $GITHUB_WORKSPACE/.github/resources/hadoop/* $HADOOP_DIR/etc/hadoop &&
  $HADOOP_DIR/bin/hdfs namenode -format &&
  $HADOOP_DIR/sbin/start-dfs.sh &&
  echo "configure_hadoop successful"
}

setup_paths() {
  echo "export JAVA_HOME=/usr/java/latest" > $HADOOP_ENV
  echo "export PATH=$HADOOP_DIR/bin:$PATH" >> $HADOOP_ENV
  echo "export LD_LIBRARY_PATH=$HADOOP_DIR/lib:$LD_LIBRARY_PATH" >> $HADOOP_ENV
  HADOOP_CP=`$HADOOP_DIR/bin/hadoop classpath --glob`
  echo "export CLASSPATH=$HADOOP_CP" >> $HADOOP_ENV
  echo "setup_paths successful"
}

install_hadoop() {
  install_prereqs
  # Check for Hadoop itself, not $HADOOP_ENV, as only $HADOOP_DIR is restored from the workflow cache
  if [[ -x $HADOOP_DIR/bin/hdfs ]]; then
    echo "Found cached Hadoop install at $HADOOP_DIR"
  else
    download_hadoop
  fi &&
    setup_paths &&
    cp -fr $GITHUB_WORKSPACE/.github/resources/hadoop/* $HADOOP_DIR/etc/hadoop &&
    mkdir -p $HADOOP_DIR/logs &&
    export HADOOP_ROOT_LOGGER=ERROR,console &&
    echo "install_hadoop successful"
  if [[ $? != 0 ]]; then
    echo "Hadoop did not install successfully. Aborting"
    exit 1
  fi
  source $HADOOP_ENV &&
    configure_hadoop &&
    echo "Install Hadoop SUCCESSFUL"
}

echo "INSTALL_DIR=$INSTALL_DIR"
install_hadoop
