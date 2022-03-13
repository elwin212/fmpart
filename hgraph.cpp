#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
// #include <process.h>
#include "hgraph.h"
#include <bits/stdc++.h>

using namespace std;

long TT = (long)time(NULL);
extern bool MODE;
extern int initcut;
//#define debug1

//#define debug
//#define bb

void error(char* s, char* s2="")
{
   cout << s <<' '<<s2<<'\n';
   exit(1);
}

cell::cell(int num,int side)
{
   number=num;
   block=side;
   gainbucket=NULL;
   gain=0;
   gainfromlock=LOCKED;
}

net::net(int num)
{
   number=num;
   unlock[0]=unlock[1]=0;
   lock[0]=lock[1]=0;
}


void parthgraph::getgraph(char *filename)
{
   int side=0, count, padoffset, temp;
   node *nnode;
   char type, state, line[80];
   FILE *fr;

   //srand(SEED);
   srand(TT);

   balance[0]=0;
   balance[1]=0;

   fr = fopen(filename, "r");
   if(!fr){
    error("cannot open input file ",filename);
   }
   //cout << filename << endl;

   fscanf(fr, "%*d %*d %d %d %d", &numnets, &numcells, &padoffset); //for scanf function %*d allow us to ignore

   //cout << numcells << " " << numnets << " " << padoffset << " " << endl;
   cells=(cell*)new cell[numcells];
   if(!cells) {
      error("memory allocation error");
   }

   minsize[0] = minsize[1] = int(ceil(numcells*RATIO));
   //cout << numcells << "  "<< minsize[0] << endl;

   int size = numcells >> 1;  //shift to right for one bit
   //cout << size << endl;
   int max[2] = {size, numcells - size};  // max size of each part
   // ***read in nodes
   for(count=0;count<numcells;count++)
   {
	  cells[count].number=count;
	  side=(rand() & 256)?1:0;
	  if(balance[side]==max[side]){
        side=other(side); // check whether both side is nearly equal
	  }
	  cells[count].block=side;
	  cells[count].gainbucket=NULL;
	  cells[count].gain=0;
	  cells[count].gainfromlock=LOCKED;
	  balance[side]++;
   }
   // ***read in nets
   nets=(net*)new net[numnets];
   if(!nets){
     error("memory allocation error");
   }
   count = -1;
   while (fgets(line, 80, fr) != NULL)
   {
	  if (sscanf(line, "%c%d %c", &type, &temp, &state) != 3){  // skip the first 5 line (# of cells, # of nets...
         continue;
	  }

	  if (type == 'p'){
         temp += padoffset;
	  }

	  if (state == 's') {
		 count++;
		 nets[count].number=0;
		 nets[count].lock[0]=0;
		 nets[count].lock[1]=0;
		 nets[count].unlock[0]=0;
		 nets[count].unlock[1]=0;
	  }
	  nnode=new node(count);
	  if(!nnode) error("memory allocation error");
	  cells[temp].first.addhead(nnode);
//	  cells[temp].first.addtail(nnode);
	  nnode=new node(temp);
	  if(!nnode) error("memory allocation error");
	  nets[count].first.addhead(nnode);
//	  nets[count].first.addtail(nnode);
	  nets[count].unlock[cells[temp].block]++;
	  nets[count].number++;
   }
   pmax=0;
   for(count=0;count<numcells;count++)
   {
	  if(cells[count].first.length>pmax)
		 pmax=cells[count].first.length;
   }
   buckets[0]=(ll) new LL[2*pmax+1];
   buckets[1]=(ll) new LL[2*pmax+1];
   if(!buckets[0] || !buckets[1]) error("memory allocation error");
   #ifdef bb
   printNodes();
   #endif
   initcut = cutset();
   if(MODE == false){
    printNodes();
   }
   cout << "Init Cutset: " << initcut << endl;
}


long parthgraph::currTime(){
    return TT;
}

void parthgraph::printNodes(){
    cout << "\tBisection1" << "\t  " << "\tBisection2" <<endl;
    for(int i =0; i < numcells; ++i){
        if(cells[i].block==0){
            cout << "\t" << "\t\t" << "\tCell " << i << endl;
        }
        else{
            cout << "\tCell " << i <<"\t" << "\t " << endl;
        }
    }

}


void parthgraph::initgains()
{

   int cnt;
   int from,to;

   // ****calc gain of unlocked cells
   for (cnt=0;cnt<numcells;cnt++)
   {
	  if(cells[cnt].gainfromlock==LOCKED)
	  {
		 cells[cnt].first.reset();
		 from=cells[cnt].block;
		 to=other(from);
		 while(cells[cnt].first.current!=NULL)
		 {
			int nnum=cells[cnt].first.current->number;
			if( (nets[nnum].unlock[from]+nets[nnum].lock[from])==1)
			   cells[cnt].gain++;
			if( (nets[nnum].unlock[to]+nets[nnum].lock[to])==0)
			   cells[cnt].gain--;
			++cells[cnt].first;
		 }
	  }
   }

   maxgain[0]=-pmax;
   maxgain[1]=-pmax;

   node *nnode;
   for(cnt=0;cnt<numcells;cnt++)
   {
	  if(cells[cnt].gainbucket==NULL)
	  {
		 nnode=(node*)new node(cnt);
		 if(!nnode) error("memory allocation error");
		 cells[cnt].gainbucket=nnode;
	  }
	  else
		 nnode=cells[cnt].gainbucket;
	  buckets[cells[cnt].block][cells[cnt].gain+pmax].addhead(nnode);
	  if(cells[cnt].gain>maxgain[cells[cnt].block])
		 maxgain[cells[cnt].block]=cells[cnt].gain;
   }
//   if(MODE == false){
//        for(int i=0;i<numcells;++i){
//        cout << "Cell:" <<cells[i].number << "  Gain:" << cells[i].gain << "  ";
//       }
//       cout << endl;
//   }

#ifdef debug1
   printCells();
   //exit(0);
#endif
}


void parthgraph::part(void)
{
   int prefix=0; // holds prefix sum of gains
   int done=0;

   passes=0;   // num of passes

   while(!done)
   {

	  passes++;
	  swapall();

	  getprefix(prefix);

	  //   if(prefix % 2)
	  //     prefix--;

	  if(prefix<=0)
		 done=1;
	  reinit(prefix);
   }
}

void parthgraph::swapall(void)
{
   int from,to;
   cell *bestcell;
   net *curnet;
   int done=0;

   bestcell=gethighest();

   while(!done)
   {
#ifdef debug
	  printBucket();
	  for(int i=0;i<numcells;++i){
        cout << cells[i].gain << " " ;
       }
       cout << endl;
	  printf("best cell=%d(%d) gain=%d\n", bestcell->number, bestcell->block,
			 bestcell->gain);
	  char go;
	  scanf(" %c", &go);
#endif
//      printf("Moving Cell No.%d(%d) gain=%d\n", bestcell->number, bestcell->block,bestcell->gain);
//	  char go;
//	  scanf(" %c", &go);
	  from=bestcell->block;
	  to=other(from);
	  bestcell->block=to;
	  bestcell->gainfromlock=bestcell->gain;
	  //bestcell->gainbucket=NULL;
	  bestcell->first.reset(); //bestcell point to head of list
	  balance[from]--;
	  balance[to]++;
	  while(bestcell->first.current != NULL)
	  {
		 curnet=&nets[bestcell->first.current->number];

		 if(curnet->lock[to]==0)
		 {
			if(curnet->unlock[to]==0)
			   fixgain(curnet->first,1);
			else if(curnet->unlock[to]==1)
			   fixgain(curnet->first,-1,1,to);
		 }

		 curnet->unlock[from]--;
		 curnet->lock[to]++;

		 if(curnet->lock[from]==0)
		 {
			if(curnet->unlock[from]==0)
			   fixgain(curnet->first,-1);
			else if(curnet->unlock[from]==1)
			   fixgain(curnet->first,1,1,from);
		 }
		 ++bestcell->first;
//		 for(int i=0;i<numcells;++i){
//            cout << "Cell:" <<cells[i].number << "  Gain:" << cells[i].gain << "  ";
//           }
//           cout << endl;
        }
	  bestcell=gethighest();
	  if(bestcell==NULL)
		 done=1;
   }

}

void parthgraph::fixgain(LL &first,int operation,int single,int sameside)
{
   cell *cellptr;
   int block;

   first.reset();
   if(operation==1)
   {
	  while(first.current!=NULL)
	  {
		 cellptr=&cells[first.current->number];
		 if(cellptr->gainfromlock==LOCKED)
		 {
			block=cellptr->block;
			if( !single || (block==sameside) )
			{
			   if( cellptr->gain==maxgain[block] )
				  ++maxgain[block];
			   buckets[block][cellptr->gain+pmax].removenode(cellptr->gainbucket);
			   cellptr->gain++;
			   buckets[block][cellptr->gain+pmax].addhead(cellptr->gainbucket);
			}
		 }
		 ++first;
	  }
   }
   else
   {
	  while(first.current!=NULL)
	  {
		 cellptr=&cells[first.current->number];
		 if(cellptr->gainfromlock==LOCKED)
		 {
			block=cellptr->block;
			if( !single || (block==sameside) )
			{
			   if( (cellptr->gain==maxgain[block]) && (buckets[block][cellptr->gain+pmax].length==1) )
				  maxgain[block]--;
			   buckets[block][cellptr->gain+pmax].removenode(cellptr->gainbucket);
			   cellptr->gain--;
			   buckets[block][cellptr->gain+pmax].addhead(cellptr->gainbucket);
			}
		 }
		 ++first;
	  }
   }
}

cell* parthgraph::gethighest(void)
{
   cell *tmpcell=NULL;
   node *tmpnode=NULL;
   //static int swapside=0;  //*******force toggle swap
   //cout << maxgain[0] << "  " << maxgain[1] << endl;
   if( (maxgain[0]>=maxgain[1]) )
   {
	  if(balance[0]>=minsize[0])
		 //  if(!swapside)  //******force toggle swap
	  {
		 tmpnode=buckets[0][maxgain[0]+pmax].removehead();
		 if(tmpnode!=NULL)
		 {
			tmpcell=&cells[tmpnode->number];
			free.addtail(tmpnode);
		 }
		 while( (!buckets[0][maxgain[0]+pmax].length) && (maxgain[0]>-pmax) )
			maxgain[0]--;
	  }
   }
   if(tmpcell==NULL)
   {
	  if(balance[1]>=minsize[1])
		 //  if(swapside)  //*******force toggle swap
	  {
		 tmpnode=buckets[1][maxgain[1]+pmax].removehead();
		 if(tmpnode!=NULL)
		 {
			tmpcell=&cells[tmpnode->number];
			free.addtail(tmpnode);
		 }
		 while( (!buckets[1][maxgain[1]+pmax].length) && (maxgain[1]>-pmax) )
			maxgain[1]--;
	  }
		 else
		 {
			tmpnode=buckets[0][maxgain[0]+pmax].removehead();
			if(tmpnode!=NULL)
			{
			   tmpcell=&cells[tmpnode->number];
			   free.addtail(tmpnode);
			}
			while( (!buckets[0][maxgain[0]+pmax].length) && (maxgain[0]>-pmax) )
			   maxgain[0]--;
		 }
   }
   //swapside=other(swapside);  //******force toggle swap
   return tmpcell;
}

void parthgraph::getprefix(int& prenum)
{
   int max=-pmax,current=0;
   int curnum=0;

   free.reset();
   while(free.current!=NULL)
   {
	  current += cells[free.current->number].gainfromlock;
	  curnum++;
	  if(max<current)
	  {
		 max=current;
		 prenum=curnum;
	  }
	  ++free;
   }
   if(max<=0)
	  prenum=0;
}

void parthgraph::reinit(int prenum)
{
   int cnt=0;
   cell *cellptr;

   if(free.length!=numcells)
   {
	  for(cnt=0;cnt<numcells;cnt++)
	  {
		 if(cells[cnt].gainfromlock==LOCKED)
		 {
			buckets[cells[cnt].block][cells[cnt].gain+pmax].removenode(cells[cnt].gainbucket);
			cells[cnt].gain = 0;  // reset to 0 --WD
		 }
	  }
   }
   free.reset();
   while(free.current!=NULL)
   {
	  cnt++;
	  cellptr=&cells[free.current->number];
	  if(cnt>prenum)
	  {
		 cellptr->block=other(cellptr->block);
		 balance[cellptr->block]++;
		 balance[other(cellptr->block)]--;
	  }
	  cellptr->gainfromlock=LOCKED;
	  cellptr->gain=0;
	  free.removenode(cellptr->gainbucket);
	  ++free;
   }

   for(cnt=0;cnt<numnets;cnt++)
   {
	  nets[cnt].lock[0]=0;
	  nets[cnt].lock[1]=0;
	  nets[cnt].unlock[0]=0;
	  nets[cnt].unlock[1]=0;
	  // ***build cell list for each net and initialize net
	  nets[cnt].first.reset();
	  while(nets[cnt].first.current!=NULL)
	  {
		 nets[cnt].unlock[ cells[nets[cnt].first.current->number].block ]++;
		 ++nets[cnt].first;
	  }
   }

   initgains();
}

int parthgraph::cutset(void)   //define what does cutset mean, if two nodes connect to each other and in different bisection, cutset count
{
   int cutset=0;
   int count;
   for(count=0;count<numnets;count++)
   {
	  //  if( ( (nets[count].unlock[0]>0) || (nets[count].lock[0]>0) ) &&
	  //      ( (nets[count].unlock[1]>0) || (nets[count].lock[1]>0) ) )
	  nets[count].first.reset();
	  int sidezero=0;
	  int sideone=0;
	  while(nets[count].first.current!=NULL)
	  {
		 if(cells[nets[count].first.current->number].block==0)
			sidezero=1;
		 else
			sideone=1;
		 ++nets[count].first;
	  }
	  if( (sideone==1) && (sidezero==1) )
		 cutset++;
   }
   return cutset;
}

void parthgraph::print(char *str)
{
   sprintf(str,"LeftSide:%i, RightSide:%i,Passes: %i",balance[0],balance[1], passes);
}

void parthgraph::printCells()
{
   cout << "Node\tPart\tGain" << endl;
   for(int i=0; i<numcells; i++)
   {
      cout << cells[i].number <<"\t"<< cells[i].block << "\t" << cells[i].gain << endl;
   }
}

void parthgraph::printBucket()
{
   for(int part=0; part<2; part++)
   {
      printf("****** Top Bucket in part %d *****\n", part);
      buckets[part][maxgain[part]+pmax].reset();
      while(buckets[part][maxgain[part]+pmax].current!=NULL)
      {
         cout << buckets[part][maxgain[part]+pmax].current->number << "\t";
		 ++buckets[part][maxgain[part]+pmax];
      }
      cout << endl;
   }
}

