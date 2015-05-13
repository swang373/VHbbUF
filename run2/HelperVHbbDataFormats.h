typedef struct {
    int HiggsFlag;
    float mass;
    float pt;
    float eta;
    float phi;
    float dR;
    float dPhi;
    float dEta;
} HiggsInfo;

typedef struct {
    int FatHiggsFlag;
    float mass;
    float pt;
    float eta;
    float phi;
    float filteredmass;
    float filteredpt;
    float filteredeta;
    float filteredphi;
    //float dR;
    //float dPhi;
    //float dEta;
} FatHiggsInfo;

typedef struct {
    float mass;  //MT in case of W
    float pt;
    float eta;
    float phi;
} VInfo;

typedef struct {
    int run;
    int lumi;
    int event;
    int json;
} EventInfo;

typedef struct {
    float et;
    float sumet;
    float sig;
    float phi;
} METInfo;

typedef struct {
  float mht;
  float ht;  
  float sig;
  float phi;
} MHTInfo;

typedef struct {
    float mass;
    float pt;
    float eta;
    float phi;
    float status;
    float charge;
    float momid;
} genParticleInfo;

typedef struct 
{
  float mass;
  float pt;
  float wMass;
} TopInfo;
