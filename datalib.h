class cl_data{
  public:
	short ant_N;
	float format;
	short memLength;//1000
	short IntpDist;
	short trig_Add;
	float trig_level;
	short v_Unit;
	//short v_Unit_div;
	//short v_Unit_ext;
	float v_Unit_div;
	float v_Unit_ext;
	//short probe_Type;
	float probe;
	float v_Scale;
	float v_Position;
	float h_scale;
	float h_Position;
	float sampl_Period;
	//float h_scale_Old;
	//float h_Position_Old;
	float timeStamp[3];
	//char timeStamp[72];
	char date[19];
	double globaltime;
	char WaveForm[2000];// (this value depend on the Memory Length  line 3)  
};
