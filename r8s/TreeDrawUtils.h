double 		calcMaxToTip(NODETYPE* root);

int 		maxorder(NODETYPE *),
			assignY2(NODETYPE *node,  double *YcurPtr,  double yinc),
			SetOneNode(NODETYPE *),
			doQWriteFile(struct MyRecType *DataPtr)
			;

NODETYPE	
		 	*TraverseScanPt(Point localPt, int radius, NODETYPE *node),
		 	*QTreeNodes(Point localPt, int radius, NODETYPE *node),
		 	*ScanInternNodes(Point localPt, int radius, NODETYPE *node),
		 	*SearchTreeNodes(Point globalPt, int radius, NODETYPE *root);

void		Tprint(NODETYPE *),
			MyDrawString(char *s, int flag),
			DrawIntName(NODETYPE *node,Rect *contRectPtr),
			DrawHigherName(struct MyRecType *DataPtr, Rect *locContentRectPtr),
			SetInternalNodeCompact(NODETYPE *node,int mode),
			ClearCompactNodes(NODETYPE *node),
			QSetAllNodes(NODETYPE *node,int flag),
			QToggleClade(NODETYPE *node),
			doQToggleTaxon(NODETYPE *node),
			assignX(NODETYPE *,  int,  int,int width,int treeMode),
			assignY(NODETYPE *,  int,  int),
			Assign_XY_Tree(NODETYPE *root,  Rect *TreeRectPtr, int treeMode),
			MagnifyTree(NODETYPE *node, int xLeft, int yLeft, double xFactor, double yFactor),
			MaxTaxLength(NODETYPE *node),
	 		MacDrawTree(NODETYPE *node),
			DrawTree(WindowPtr w),

			SetCompactNodes(NODETYPE *),
			MakeTreeStruct(char *TreeDescriptionPtr,struct MyRecType *theDataPtr),

			Dispose(struct MyRecType *theDataPtr);

void drawRuler(Rect TreeRect,int width, int treeMode, double maxLength);
