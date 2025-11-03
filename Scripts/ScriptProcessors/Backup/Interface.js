Content.makeFrontInterface(1200, 760);


include("Neural1.js");
const var namModel1 = Engine.createNeuralNetwork("NN_1");
namModel1.loadNAMModel(Model1);

include("Neural2.js");
const namModel2 = Engine.createNeuralNetwork("NN_2");
namModel2.loadNAMModel(Model2); 



Synth.deferCallbacks(true);

// Keyswitches & muters
const var keySwitches = [67, 69]; // G3, A3
const var midiMuters = [];
reg i;

for (i = 0; i < keySwitches.length; i++)
    midiMuters[i] = Synth.getMidiProcessor("MidiMuter" + i);

// Buttons
const var artButtons = [];
const var btnNames = ["DownPick", "Alternate Pick"];

// Helpers
inline function updateButtons(index)
{
    for (i = 0; i < artButtons.length; i++)
        artButtons[i].setValue(i == index);
}

inline function changeArticulation(index)
{
    for (i = 0; i < midiMuters.length; i++)
        midiMuters[i].setAttribute(0, i != index); // mute others
    updateButtons(index);
}

// Two explicit callbacks (no component properties needed)
inline function Artic0(c, v) { if (v) changeArticulation(0); }
inline function Artic1(c, v) { if (v) changeArticulation(1); }

// Create buttons & assign the right callback (single creation)
for (i = 0; i < btnNames.length; i++)
{
    artButtons[i] = Content.addButton("ArtBtn" + i, 20 + (i * 150), 40);
    artButtons[i].set("text", btnNames[i]);
    artButtons[i].set("saveInPreset", false);
    artButtons[i].setControlCallback(i == 0 ? Artic0 : Artic1);
}

// Set fixed positions on the already-created buttons (no re-adding)
artButtons[0].set("x", 199);
artButtons[0].set("y", 137);

artButtons[1].set("x", 793);
artButtons[1].set("y", 155);

// Init selection
changeArticulation(0);
Synth.deferCallbacks(false);

// ======= MIDI =======
reg tmpNote, tmpIdx;

function onNoteOn()
{
    tmpNote = Message.getNoteNumber();
    tmpIdx  = keySwitches.indexOf(tmpNote);

    if (tmpIdx != -1)
    {
        changeArticulation(tmpIdx); // sync buttons and muters
        Message.ignoreEvent(true);  // swallow keyswitch note
    }
}

function onNoteOff()
{
    tmpNote = Message.getNoteNumber();
    if (keySwitches.indexOf(tmpNote) != -1)
        Message.ignoreEvent(true);
}

function onController() {}
function onAftertouch() {}
function onPitchWheel() {}


function onNoteOn()
{

    if (keySwitches.indexOf(Message.getNoteNumber()) != -1) //If a keyswitch was pressed
    {
        changeArticulation(keySwitches.indexOf(Message.getNoteNumber())); //Call the changeArticulation function
    }    
}
function onNoteOff()
{
	
}
 function onController()
{
	
}
 function onTimer()
{
	
}
 function onControl(number, value)
{
	
}
 