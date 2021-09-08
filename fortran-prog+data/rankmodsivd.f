c
c compute N(r,t,n) for r=0,1,..(t-1)/2 for given t and n
c
	implicit integer(a-y)
	dimension p(0:16000000)
        dimension nr(0:100)
        dimension ncr(0:100)
c
c p(j) = number of partitions of j (j=0,1,...n)
c nr(r) = N(r,t,n) for given n  (r=0,1,..(t-1))
c ncr(r) = M(r,t,n) for given n  (r=0,1,..(t-1))
c
        open(unit=1,file='inputdata2')
        open(unit=2,file='ranksave')
        open(unit=3,file='mdim')
        open(unit=4,file='ranksave2')
        open(unit=12,file='cranksave')
        open(unit=13,file='ptnsave')
        open(unit=14,file='cranksave2')
 
c
c file inputdata contains two integers; eg.
c 11
c 172904
c
	p(0)=1
	read(1,*) t
	read(1,*) ln
        read(1,*) smod
        read (1,*) sres
	write(*,*) 't = ',t
	write(*,*) 'n = ',ln
	write(*,*) 'smod = ',smod
	write(*,*) 'sres = ',sres
        t1=(t-1)/2
c-----------------------------------------------------------------
c compute p(j) mod t  for j=0,1,..,n -----------------------------
c-----------------------------------------------------------------
        cs=0
        cns=0
	do 25 i=1,ln
	p(i)=0
25	continue
        do 24 j=0,100
        nr(j)=0
24      continue
        write(*,*) "init done"
	do 100 n=1,ln
         nn1=modp(n,10000)
         if (nn1.eq.0) then
           write(*,*) "computing p(",n,") mod ",t
         end if
         zn=n
         zk1=(sqrt(24*zn+1)+1)/6
         kz=zk1
	 do 50 k=1,kz
	 nk1=n-k*(3*k+1)/2
	 nk2=n-k*(3*k-1)/2
	 if (nk1.ge.0) then
	 p(n)=p(n) - (-1)**k*p(nk1)
	 xp=modp(p(n),t)
	 p(n)=xp
	 end if
	 if (nk2.ge.0) then
	 p(n)=p(n) - (-1)**k*p(nk2)
	 xp=modp(p(n),t)
	 p(n)=xp
	 end if
50	 continue
100      continue
c        ln = j*smod + sres 
	 lj=(ln-sres)/smod
	 do 1000 jj=0,lj
            jjm=modp(jj,100)
            if (jjm.eq.0) then
               write(*,*) jj,lj
	    endif
c-------------------------------------------------
c reset nr array to zero
c
           do 524 j=0,40
             nr(j)=0
524        continue
c-------------------------------------------------
c         n=ln
c         zn=ln
	 n=jj*smod+sres
	 write(13,*) p(n)
         zn=n
         zk1=(sqrt(24*zn+1)+1)/6
         kz=zk1
c         write(*,*) "kz=",kz
	 do 200 k=1,kz+1
         km=modp(k,100)
         if (km.eq.0) then
c           write(*,*) "k=",k,",  lastk=",kz
	 end if
	 tt2=n-k*(3*k+1)/2
	 tt1=n-k*(3*k-1)/2
	 if (tt1.ge.0) then
              m1=tt1/k
              do 140 m=0,m1
              xx1=(-1)**(k-1)*p(tt1-k*m)
              if (m.eq.0) then
                 nr(0)=nr(0)+xx1
                 nr(0)=modp(nr(0),t)
              else
                  mr=modp(m,t)
                  if (mr.eq.0) then
                     nr(0)=nr(0)+2*xx1
                     nr(0)=modp(nr(0),t)
                  else
                     nr(mr)=nr(mr)+xx1
                     nr(mr)=modp(nr(mr),t)
                     nr(t-mr)=nr(t-mr)+xx1
                     nr(t-mr)=modp(nr(t-mr),t)
                  end if
              end if
140         continue
	 end if
           if (tt2.ge.0) then
              m1=tt2/k
              do 150 m=0,m1
              xx1=(-1)**(k)*p(tt2-k*m)
              if (m.eq.0) then
                 nr(0)=nr(0)+xx1
                 nr(0)=modp(nr(0),t)
              else
                  mr=modp(m,t)
                  if (mr.eq.0) then
                     nr(0)=nr(0)+2*xx1
                     nr(0)=modp(nr(0),t)
                  else
                     nr(mr)=nr(mr)+xx1
                     nr(mr)=modp(nr(mr),t)
                     nr(t-mr)=nr(t-mr)+xx1
                     nr(t-mr)=modp(nr(t-mr),t)
                  end if
              end if
150         continue
         end if
200      continue
c---------------------------------------------------------
c begin crank loop
c-------------------------------------------------
c reset ncr array to zero
c
           do 1524 j=0,40
             ncr(j)=0
1524        continue
c-------------------------------------------------
c         n=ln
c         zn=ln
	 n=jj*smod+sres
         zn=n
         zk1=(sqrt(8*zn+1)+1)/2
         kz=zk1
c         write(*,*) "kz=",kz
	 do 1200 k=1,kz+1
         km=modp(k,100)
         if (km.eq.0) then
c           write(*,*) "k=",k,",  lastk=",kz
	 end if
	 tt2=n-k*(k+1)/2
	 tt1=n-k*(k-1)/2
	 if (tt1.ge.0) then
              m1=tt1/k
              do 1140 m=0,m1
              xx1=(-1)**(k-1)*p(tt1-k*m)
              if (m.eq.0) then
                 ncr(0)=ncr(0)+xx1
                 ncr(0)=modp(ncr(0),t)
              else
                  mr=modp(m,t)
                  if (mr.eq.0) then
                     ncr(0)=ncr(0)+2*xx1
                     ncr(0)=modp(ncr(0),t)
                  else
                     ncr(mr)=ncr(mr)+xx1
                     ncr(mr)=modp(ncr(mr),t)
                     ncr(t-mr)=ncr(t-mr)+xx1
                     ncr(t-mr)=modp(ncr(t-mr),t)
                  end if
              end if
1140         continue
	 end if
           if (tt2.ge.0) then
              m1=tt2/k
              do 1150 m=0,m1
              xx1=(-1)**(k)*p(tt2-k*m)
              if (m.eq.0) then
                 ncr(0)=ncr(0)+xx1
                 ncr(0)=modp(ncr(0),t)
              else
                  mr=modp(m,t)
                  if (mr.eq.0) then
                     ncr(0)=ncr(0)+2*xx1
                     ncr(0)=modp(ncr(0),t)
                  else
                     ncr(mr)=ncr(mr)+xx1
                     ncr(mr)=modp(ncr(mr),t)
                     ncr(t-mr)=ncr(t-mr)+xx1
                     ncr(t-mr)=modp(ncr(t-mr),t)
                  end if
              end if
1150         continue
         end if
1200	continue
c---------------------------------------------------------
c         do 900 a=0,t1-1
c          ff=ff+(nr(a)-nr(a+1))**2
c900	continue
c         if (ff.eq.0) then
c           cs=cs+1
c         else
c           cns=cns+1
c         endif
cc         write(*,*) "N(r,",t,",",n,"), r=0,1,..",t1,"cs cns jjm:"
cc         write(*,*) (nr(j),j=0,t1),"   ",cs,cns,modp(jj,t)
         cd=1
         do 950 j=1,t1
           dd=nr(0)-nr(j)
           if (dd.ne.0) then
              cd=0
           endif
950	 continue
         write(2,*) (nr(j),j=0,t1)
         if (cd.eq.1) then
           write(4,*) n,(n-sres)/smod,nr(0)
         endif
         ff=0
         xcd=1
         do 1950 j=1,t1
           xdd=ncr(0)-ncr(j)
           if (xdd.ne.0) then
              xcd=0
           endif
1950	 continue
         write(12,*) (ncr(j),j=0,t1)
         if (xcd.eq.1) then
           write(14,*) n,(n-sres)/smod,ncr(0)
         endif
         ff=0
1000	continue
        write(3,*) t      
        write(3,*) t1
        write(3,*) smod
        write(3,*) sres
        write(3,*) lj+1
	stop
	end
	function modp(n,t)
	implicit integer(a-z)
	if (n.ge.0) then
	modp=mod(n,t)
	else
	x1=mod(-n,t)
	if (x1.eq.0) then
	modp=x1
	else
	modp=t-x1
	end if
	end if
	return
	end
