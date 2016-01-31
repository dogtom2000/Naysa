



var Earth = {
    gravPara: 398600,       //km3/s2
    radius: 6371,           //km
    pressure: 1,            //atm
    atmScale: 7,            //dimensionless
    atmHeight: 140,         //km
    atmWeight: 28.97,       //g/mol
    velEquator: 465.1,      //m/s
}

var Rocket = {
    massShip: 1000,         //kg
    massFuel: 20000,        //kg
    thrust: 250000,         //N
    Isp: 1000,              //s
    cD: 0.2,                //dimensionless
    diameter: 5,            //m
    velocity: [0,0],        //m/s [x,z]
    position: [0,0]         //km [x,z]
}


var orbit = 180;            //km what you want your orbit to be

var gravBase = 9.81;        //m/s2
var gasConst = 8.314;       //m3Pa/K/mol
var rocketArea = Math.PI / 4 * Math.pow(Rocket.diameter,2); //m2
var stageFuel = Rocket.massFuel; //kg
var apoFlag = 0;
var apoapsis = 0;



function burn(deltaTime){
    
    var vel_X_old = Rocket.velocity[0];
    var vel_Z_old = Rocket.velocity[1];
    var pos_X_old = Rocket.position[0];
    var pos_Z_old = Rocket.position[1];
    var massFuel_old = Rocket.massFuel;
    var apo_old = apoapsis;
    
    
    var velocityMag = Math.pow((Math.pow(Rocket.velocity[0],2) + Math.pow(Rocket.velocity[1],2)),0.5);                    //m/s

    if (Rocket.position[1] < Earth.atmHeight){
        var atmDens = 0.042295 * Earth.pressure * Math.exp(-Rocket.position[1] / Earth.atmScale) * Earth.atmWeight;            //kg/m3 (Temp = 288.15 K)
        var aDrag = Rocket.cD * atmDens * Math.pow(velocityMag, 2) * rocketArea / 2 / (Rocket.massShip + Rocket.massFuel);      //m/s2
    } else {
        var aDrag = 0;      //m/s2
    }
    
    if (velocityMag > 0){
    var aDrag_Z = aDrag * Rocket.velocity[1]/velocityMag;                                                                  //m/s2
    var aDrag_X = aDrag * Rocket.velocity[0]/velocityMag;                                                                  //m/s2   

    } else {
        var aDrag_Z = 0;
        var aDrag_X = 0;
    }

    var aCent_Z = Math.pow(Rocket.velocity[0] + Earth.velEquator,2) / (Earth.radius + Rocket.position[1]) / 1000;                            //m/s2
    var aGrav_Z = Earth.gravPara / Math.pow((Earth.radius + Rocket.position[1]),2) * 1000;                                 //m/s2
    
    apoapsis = Math.pow(Rocket.velocity[1],2) / 2 / (aGrav_Z) / 1000 + Rocket.position[1];                                        //km


    if (apoFlag === 0){
        
        var pitch = Math.PI / 2;
        
    } else if (Rocket.velocity[1] <= 0){
        
       var pitch = Math.asin((aDrag_Z + aGrav_Z) / Rocket.thrust * (Rocket.massShip + Rocket.massFuel));
        
    } else {
        var pitch = Math.asin((aDrag_Z) / Rocket.thrust * (Rocket.massShip + Rocket.massFuel));
    }

    var aThrust_Z = Rocket.thrust * Math.sin(pitch) / (Rocket.massShip + Rocket.massFuel);                                  //m/s2
    var aThrust_X = Rocket.thrust * Math.cos(pitch) / (Rocket.massShip + Rocket.massFuel);                                  //m/s2
    
    Rocket.velocity[0] += (aThrust_X - aDrag_X) * deltaTime;
    Rocket.velocity[1] += (aThrust_Z - aDrag_Z + aCent_Z - aGrav_Z) * deltaTime;
    Rocket.position[1] += (vel_Z_old + Rocket.velocity[1]) / 2 * deltaTime / 1000;
    
    Rocket.massFuel -= Rocket.thrust / Rocket.Isp / gravBase * deltaTime;
    
    if (apoapsis > orbit && apoFlag === 0){
        apoFlag = 1;
        
        //console.log("here", Rocket.velocity[0], Rocket.velocity[1], vel_X_old, vel_Z_old, apoapsis, apo_old)
        
        Rocket.velocity[0] = (Rocket.velocity[0] - vel_X_old) * (orbit - apoapsis) / (apoapsis - apo_old) + vel_X_old;
        Rocket.velocity[1] = (Rocket.velocity[1] - vel_Z_old) * (orbit - apoapsis) / (apoapsis - apo_old) + vel_Z_old;
        Rocket.position[0] = (Rocket.position[0] - pos_X_old) * (orbit - apoapsis) / (apoapsis - apo_old) + pos_X_old;
        Rocket.position[1] = (Rocket.position[1] - pos_Z_old) * (orbit - apoapsis) / (apoapsis - apo_old) + pos_Z_old;
        Rocket.massFuel = (Rocket.massFuel - massFuel_old) * (orbit - apoapsis) / (apoapsis - apo_old) + massFuel_old;
        
        Rocket.position[1] = orbit - Math.pow(Rocket.velocity[1],2) / 2 / (aGrav_Z) / 1000;
        apoapsis = Math.pow(Rocket.velocity[1],2) / 2 / (aGrav_Z) / 1000 + Rocket.position[1]; 
        //console.log("here", Rocket.velocity[0], Rocket.velocity[1])
        velocityMag = Math.pow((Math.pow(Rocket.velocity[0],2) + Math.pow(Rocket.velocity[1],2)),0.5);
    }
    

    console.log(Rocket.velocity[1])
    //console.log(Math.log(Rocket.massShip + Rocket.massFuel)*9.81*Rocket.Isp, velocityMag)
    
    if (aCent_Z > aGrav_Z){
        Rocket.velocity[1] = 0;
        Rocket.velocity[0] = Math.pow((aGrav_Z * (Earth.radius + Rocket.position[1]) * 1000),0.5)
        Rocket.position[1] = orbit;
        velocityMag = Math.pow((Math.pow(Rocket.velocity[0],2) + Math.pow(Rocket.velocity[1],2)),0.5);
        //console.log(Math.log(Rocket.massShip + Rocket.massFuel)*9.81*Rocket.Isp, velocityMag)
        console.log(apoapsis)
        return 0;
    } else {
        return 1;
    }
    
}


for (var i = 0; i < 20; i++){
 if (burn(stageFuel / (Rocket.thrust / Rocket.Isp / gravBase) / 20) === 0){
     break;
 }
};


