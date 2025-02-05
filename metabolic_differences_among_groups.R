
write.table(GSE25066_patient, file = "GSE25066_patient.txt", sep = "\t",
            row.names = FALSE)


library(tidyverse)
library(ggpubr)

GSE25066_patient = t(GSE25066_patient)
GSE25066_patient = as.data.frame(GSE25066_patient)
##Order 
clusters_id_GSE25066 = c("GSM615110",
                         "GSM615384",
                         "GSM615739",
                         "GSM615147",
                         "GSM615167",
                         "GSM615127",
                         "GSM615125",
                         "GSM615172",
                         "GSM615120",
                         "GSM615114",
                         "GSM615122",
                         "GSM615158",
                         "GSM615166",
                         "GSM615170",
                         "GSM615173",
                         "GSM615162",
                         "GSM615646",
                         "GSM615168",
                         "GSM615176",
                         "GSM615142",
                         "GSM615144",
                         "GSM615163",
                         "GSM615174",
                         "GSM615306",
                         "GSM615795",
                         "GSM615129",
                         "GSM615165",
                         "GSM615150",
                         "GSM615097",
                         "GSM615133",
                         "GSM615828",
                         "GSM615134",
                         "GSM615138",
                         "GSM615177",
                         "GSM615299",
                         "GSM615314",
                         "GSM615100",
                         "GSM615820",
                         "GSM615760",
                         "GSM615758",
                         "GSM615761",
                         "GSM615768",
                         "GSM615772",
                         "GSM615781",
                         "GSM615157",
                         "GSM615131",
                         "GSM615148",
                         "GSM615159",
                         "GSM615757",
                         "GSM615742",
                         "GSM615272",
                         "GSM615389",
                         "GSM615320",
                         "GSM615373",
                         "GSM615119",
                         "GSM615287",
                         "GSM615650",
                         "GSM615399",
                         "GSM615819",
                         "GSM615815",
                         "GSM615794",
                         "GSM615804",
                         "GSM615330",
                         "GSM615382",
                         "GSM615656",
                         "GSM615223",
                         "GSM615338",
                         "GSM615640",
                         "GSM615344",
                         "GSM615385",
                         "GSM615257",
                         "GSM615233",
                         "GSM615386",
                         "GSM615318",
                         "GSM615316",
                         "GSM615752",
                         "GSM615747",
                         "GSM615268",
                         "GSM615298",
                         "GSM615188",
                         "GSM615294",
                         "GSM615280",
                         "GSM615289",
                         "GSM615256",
                         "GSM615283",
                         "GSM615404",
                         "GSM615279",
                         "GSM615789",
                         "GSM615259",
                         "GSM615180",
                         "GSM615273",
                         "GSM615254",
                         "GSM615292",
                         "GSM615269",
                         "GSM615286",
                         "GSM615281",
                         "GSM615245",
                         "GSM615253",
                         "GSM615266",
                         "GSM615247",
                         "GSM615258",
                         "GSM615293",
                         "GSM615728",
                         "GSM615737",
                         "GSM615666",
                         "GSM615680",
                         "GSM615369",
                         "GSM615376",
                         "GSM615357",
                         "GSM615368",
                         "GSM615371",
                         "GSM615317",
                         "GSM615354",
                         "GSM615202",
                         "GSM615397",
                         "GSM615307",
                         "GSM615326",
                         "GSM615321",
                         "GSM615672",
                         "GSM615213",
                         "GSM615211",
                         "GSM615212",
                         "GSM615181",
                         "GSM615276",
                         "GSM615649",
                         "GSM615632",
                         "GSM615677",
                         "GSM615676",
                         "GSM615393",
                         "GSM615716",
                         "GSM615290",
                         "GSM615644",
                         "GSM615336",
                         "GSM615324",
                         "GSM615641",
                         "GSM615651",
                         "GSM615194",
                         "GSM615661",
                         "GSM615798",
                         "GSM615639",
                         "GSM615648",
                         "GSM615209",
                         "GSM615732",
                         "GSM615260",
                         "GSM615380",
                         "GSM615771",
                         "GSM615779",
                         "GSM615334",
                         "GSM615660",
                         "GSM615112",
                         "GSM615348",
                         "GSM615191",
                         "GSM615681",
                         "GSM615310",
                         "GSM615347",
                         "GSM615674",
                         "GSM615242",
                         "GSM615787",
                         "GSM615248",
                         "GSM615741",
                         "GSM615251",
                         "GSM615715",
                         "GSM615764",
                         "GSM615264",
                         "GSM615288",
                         "GSM615184",
                         "GSM615193",
                         "GSM615189",
                         "GSM615250",
                         "GSM615135",
                         "GSM615182",
                         "GSM615197",
                         "GSM615215",
                         "GSM615139",
                         "GSM615754",
                         "GSM615664",
                         "GSM615271",
                         "GSM615274",
                         "GSM615265",
                         "GSM615355",
                         "GSM615670",
                         "GSM615684",
                         "GSM615706",
                         "GSM615689",
                         "GSM615702",
                         "GSM615712",
                         "GSM615687",
                         "GSM615694",
                         "GSM615699",
                         "GSM615240",
                         "GSM615244",
                         "GSM615691",
                         "GSM615696",
                         "GSM615746",
                         "GSM615668",
                         "GSM615671",
                         "GSM615300",
                         "GSM615363",
                         "GSM615714",
                         "GSM615190",
                         "GSM615667",
                         "GSM615638",
                         "GSM615185",
                         "GSM615657",
                         "GSM615727",
                         "GSM615733",
                         "GSM615394",
                         "GSM615372",
                         "GSM615375",
                         "GSM615378",
                         "GSM615390",
                         "GSM615312",
                         "GSM615387",
                         "GSM615637",
                         "GSM615145",
                         "GSM615102",
                         "GSM615152",
                         "GSM615136",
                         "GSM615149",
                         "GSM615106",
                         "GSM615151",
                         "GSM615118",
                         "GSM615762",
                         "GSM615763",
                         "GSM615156",
                         "GSM615346",
                         "GSM615098",
                         "GSM615146",
                         "GSM615099",
                         "GSM615164",
                         "GSM615116",
                         "GSM615117",
                         "GSM615160",
                         "GSM615178",
                         "GSM615776",
                         "GSM615325",
                         "GSM615396",
                         "GSM615104",
                         "GSM615113",
                         "GSM615143",
                         "GSM615140",
                         "GSM615107",
                         "GSM615109",
                         "GSM615141",
                         "GSM615748",
                         "GSM615335",
                         "GSM615643",
                         "GSM615364",
                         "GSM615291",
                         "GSM615137",
                         "GSM615333",
                         "GSM615154",
                         "GSM615175",
                         "GSM615327",
                         "GSM615773",
                         "GSM615169",
                         "GSM615171",
                         "GSM615686",
                         "GSM615700",
                         "GSM615711",
                         "GSM615713",
                         "GSM615708",
                         "GSM615743",
                         "GSM615750",
                         "GSM615688",
                         "GSM615749",
                         "GSM615685",
                         "GSM615697",
                         "GSM615654",
                         "GSM615751",
                         "GSM615635",
                         "GSM615308",
                         "GSM615662",
                         "GSM615811",
                         "GSM615340",
                         "GSM615343",
                         "GSM615631",
                         "GSM615717",
                         "GSM615339",
                         "GSM615360",
                         "GSM615297",
                         "GSM615658",
                         "GSM615319",
                         "GSM615361",
                         "GSM615805",
                         "GSM615818",
                         "GSM615653",
                         "GSM615673",
                         "GSM615231",
                         "GSM615203",
                         "GSM615337",
                         "GSM615285",
                         "GSM615208",
                         "GSM615187",
                         "GSM615277",
                         "GSM615659",
                         "GSM615183",
                         "GSM615759",
                         "GSM615342",
                         "GSM615755",
                         "GSM615722",
                         "GSM615374",
                         "GSM615383",
                         "GSM615352",
                         "GSM615204",
                         "GSM615278",
                         "GSM615719",
                         "GSM615740",
                         "GSM615331",
                         "GSM615725",
                         "GSM615721",
                         "GSM615388",
                         "GSM615724",
                         "GSM615123",
                         "GSM615153",
                         "GSM615669",
                         "GSM615332",
                         "GSM615356",
                         "GSM615636",
                         "GSM615807",
                         "GSM615777",
                         "GSM615823",
                         "GSM615826",
                         "GSM615132",
                         "GSM615775",
                         "GSM615652",
                         "GSM615825",
                         "GSM615735",
                         "GSM615731",
                         "GSM615249",
                         "GSM615395",
                         "GSM615305",
                         "GSM615821",
                         "GSM615350",
                         "GSM615353",
                         "GSM615196",
                         "GSM615678",
                         "GSM615367",
                         "GSM615778",
                         "GSM615206",
                         "GSM615679",
                         "GSM615311",
                         "GSM615709",
                         "GSM615665",
                         "GSM615790",
                         "GSM615675",
                         "GSM615302",
                         "GSM615663",
                         "GSM615304",
                         "GSM615718",
                         "GSM615736",
                         "GSM615267",
                         "GSM615201",
                         "GSM615200",
                         "GSM615216",
                         "GSM615205",
                         "GSM615227",
                         "GSM615108",
                         "GSM615121",
                         "GSM615634",
                         "GSM615642",
                         "GSM615633",
                         "GSM615645",
                         "GSM615261",
                         "GSM615246",
                         "GSM615366",
                         "GSM615792",
                         "GSM615655",
                         "GSM615802",
                         "GSM615214",
                         "GSM615161",
                         "GSM615341",
                         "GSM615365",
                         "GSM615096",
                         "GSM615370",
                         "GSM615232",
                         "GSM615243",
                         "GSM615683",
                         "GSM615814",
                         "GSM615345",
                         "GSM615329",
                         "GSM615351",
                         "GSM615228",
                         "GSM615179",
                         "GSM615210",
                         "GSM615199",
                         "GSM615237",
                         "GSM615362",
                         "GSM615207",
                         "GSM615349",
                         "GSM615103",
                         "GSM615217",
                         "GSM615682",
                         "GSM615186",
                         "GSM615198",
                         "GSM615734",
                         "GSM615756",
                         "GSM615192",
                         "GSM615647",
                         "GSM615693",
                         "GSM615234",
                         "GSM615236",
                         "GSM615238",
                         "GSM615230",
                         "GSM615745",
                         "GSM615705",
                         "GSM615710",
                         "GSM615703",
                         "GSM615707",
                         "GSM615695",
                         "GSM615701",
                         "GSM615690",
                         "GSM615698",
                         "GSM615704",
                         "GSM615224",
                         "GSM615692",
                         "GSM615767",
                         "GSM615765",
                         "GSM615766",
                         "GSM615780",
                         "GSM615813",
                         "GSM615817",
                         "GSM615769",
                         "GSM615770",
                         "GSM615784",
                         "GSM615824",
                         "GSM615782",
                         "GSM615816",
                         "GSM615827",
                         "GSM615262",
                         "GSM615313",
                         "GSM615226",
                         "GSM615744"
)

num_clusters_GSE25066 = c("Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 1",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 2",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3",
                          "Grupo 3")



GSE25066_patient = GSE25066_patient[match(clusters_id_GSE25066, GSE25066_patient$patient_id),]



GSE25066_patient = cbind(num_clusters_GSE25066, GSE25066_patient)

GSE25066_patient = transform(GSE25066_patient, 
                             MYC = as.numeric(MYC))


#################################################################################################
ggboxplot(GSE25066_patient, x = "num_clusters_GSE25066",
          y = c("MYC", "MYC", "MYC"),
          combine = TRUE,
          color = "num_clusters_GSE25066", palette = "jco",
          ylab = "Expression", 
          add = "jitter",                              # Add jittered points
          add.params = list(size = 0.1, jitter = 0.2)  # Point size and the amount of jittering
)

################################################################################################


library("survminer")


GSE25066_patient = transform(GSE25066_patient, 
                             MYC = as.numeric(MYC))

# Perform the test
MYC_p_adj_GSE25066 = compare_means(MYC~ num_clusters_GSE25066,  data = GSE25066_patient,
                                                  ref.group = ".all.",
                                                  method = "t.test", 
                                                  p.adjust.method = "bonferroni")

MYC_p_adj_GSE25066


ggboxplot(GSE25066_patient, x = "num_clusters_GSE25066", y = "MYC", color = "num_clusters_GSE25066", 
          add = "jitter", 
          legend = "none", 
          palette = "jco",
          xlab = "Grupos metabólicos",
          ylab = expression(~bolditalic("MYC"))) +
  stat_compare_means(method = "anova", 
                     label.y = 1.40) +
  theme(axis.title.y = element_text(size = 19, face = "italic", 
                                    margin = margin(t = 0, r = 18, b = 0, l = 0)),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 18, r = 0, b = 0, l = 0))) +# Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all

GSE= MYC_GROUPS_GSE25066







