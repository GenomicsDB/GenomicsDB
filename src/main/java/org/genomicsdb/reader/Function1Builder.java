package org.genomicsdb.reader;

import scala.Function1;

public class Function1Builder<T,R>{

  public static interface Function1Apply<T,R> {
    R apply(T t);
  }

  private Function1Apply<T,R> lambda;
  
  private Function1<T,R> function;
  
  public Function1Builder(Function1Apply<T,R> lambda){
    this.lambda = lambda;
    this.function = new Function1<T,R>(){
      @Override
      public R apply(T t){
        return Function1Builder.this.lambda.apply(t);
      }
    };
  }

  public Function1<T, R> getFunction(){
    return this.function;
  }
}
