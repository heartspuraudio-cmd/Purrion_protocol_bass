 const var Down = Synth.getSampler("Down");
 Down.enableRoundRobin(false);
 
 reg counter = 0;
 
 const var Alt = Synth.getSampler("Alt");
 Alt.enableRoundRobin(false);
 

 

 function onNoteOn()
{

Down.setActiveGroup(Math.randInt(1, 7));

Alt.setActiveGroup(Math.randInt(1, 7));


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
 