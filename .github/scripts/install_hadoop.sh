#!/bin/bash

# Install hadoop
# Installation relies on finding JAVA_HOME@/usr/java/latest as a prerequisite

INSTALL_DIR=${INSTALL_DIR:-/usr}
USER=`whoami`

HADOOP=hadoop-${HADOOP_VER:-3.2.1}
HADOOP_DIR=${INSTALL_DIR}/$HADOOP
HADOOP_ENV=$HADOOP_DIR/hadoop.env

install_prereqs() {
  if [[ -f /usr/java/latest ]]; then
    echo "/usr/java/latest found"
      sudo rm /usr/java/latest
  fi
  if [[ ! -z $JAVA_HOME ]]; then
    sudo mkdir -p /usr/java
    sudo ln -s $JAVA_HOME /usr/java/latest
  else 
    sudo apt install openjdk-8-jre-headless
    sudo ln -s /usr/lib/jvm/java-1.8.0-openjdk-amd64/ /usr/java/latest
  fi
  echo "install_prereqs successful"
}

download_hadoop() {
  wget -q http://www-eu.apache.org/dist/hadoop/common/$HADOOP/$HADOOP.tar.gz &&
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
  ssh-keygen -t rsa -b 4096 -N '' -f ~/.ssh/id_rsa &&
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
  if [[ ! -f $HADOOP_ENV ]]; then
    download_hadoop &&
      setup_paths &&
      cp -fr $GITHUB_WORKSPACE/.github/resources/hadoop/* $HADOOP_DIR/etc/hadoop &&
      mkdir -p $HADOOP_DIR/logs &&
      export HADOOP_ROOT_LOGGER=ERROR,console
  fi
  source $HADOOP_ENV &&
    configure_hadoop &&
    echo "Install Hadoop SUCCESSFUL"
}

echo "INSTALL_DIR=$INSTALL_DIR"
echo "INSTALL_TYPE=$INSTALL_TYPE"
install_hadoop
