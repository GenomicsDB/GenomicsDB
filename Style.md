## GenomicsDB Style Guide

For Java code, roughly adhere to [Google Java Style](https://google.github.io/styleguide/javaguide.html) and for C/C++ roughly adhere to [Google C++ Style](https://google.github.io/styleguide/cppguide.html) for consistency/readabilty 

GenomicsDB Example Rules:

```
Use spaces instead of tabs.
Use 2 spaces for indenting.
Add brackets even for one line blocks e.g. 
        if (x>0)
           do_foo();
 should ideally be 
       if (x>0) {
         do_foo();
       }
Pad header e.g.
        if(x>0) should be if (x>0)
        while(x>0) should be while (x>0)
One half indent for class modifiers.
```
