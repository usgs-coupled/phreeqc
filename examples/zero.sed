#!/bin/sh
sed "s~ \-0\.00 ~  0\.00 ~g" < $1 | sed "s~ \-0\.000 ~  0\.000 ~g" | sed "s~ \-0\.0000 ~  0\.0000 ~g" | \
  sed "s~ \-0\.00	~  0\.00	~g" | sed "s~ \-0\.000	~  0\.000	~g" | sed "s~ \-0\.0000	~  0\.0000	~g" | \
  sed "s~ \-0\.00$~  0\.00~" | sed "s~ \-0\.000$~  0\.000~" | sed "s~ \-0\.0000$~  0\.0000~" | \
  sed "s~e\-00~e\+00~g" > t
mv t $1
