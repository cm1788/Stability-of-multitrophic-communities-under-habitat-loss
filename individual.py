


from configure import MAX_RESOURCE, MIN_RESOURCE, LIVING_EXPEND, SYNTHESIS_ABILITY, MUTUALISTIC_FRACTION, OMNIVORY_TRADEOFF
from configure import MATING_RESOURCE, MATING_ENERGY, EFFICIENCY_TRANSFER, HERBIVORY_EFFICIENCY, HERBIVORY_FRACTION
from configure import MUTUALISTIC_EFFICIENCY, MUTUALISTIC_COOLING, MIN_MUTUALISTIC_EFF, MUTUALISTIC_LOSS, MUT_PRODUCER_LOSS

from random import randint

class Individual():
    
    def __init__(self, species_id):
        self.resource = float(randint(2*MIN_RESOURCE, MAX_RESOURCE))
        self.mutualist = False
        self.species_id = species_id
        self.mutualistic_producer = False
    
    def become_mutualist(self):
        self.mutualist = True
        self.current_host = None
    
    def become_mutualistic_producer(self):
        self.mutualistic_producer = True
    
    def live(self, producer):
#        producer = False
        #if not producer:
        if self.mutualist:
            lost_resource = MUTUALISTIC_LOSS*LIVING_EXPEND
        elif self.mutualistic_producer:
            lost_resource = MUT_PRODUCER_LOSS*LIVING_EXPEND
        else:
            lost_resource = LIVING_EXPEND
          
        self.resource -= self.resource*lost_resource
        #else:
        #    if self.mutualist:
        #        print 'I am a mutualist not losing any resource'
        
        if self.resource < MIN_RESOURCE:
            return False
        
        return True
    
    def eat(self, prey, herbivorous=False, omnivore=False):   
        if herbivorous:
            if self.mutualist:
                lost_biomass_prey = MAX_RESOURCE*HERBIVORY_FRACTION
                #print 'Mutualist feeding... gain = ', lost_biomass_prey
            else:
                lost_biomass_prey = MAX_RESOURCE*HERBIVORY_FRACTION
            
            if prey.resource < lost_biomass_prey:
                lost_biomass_prey = prey.resource
                     
                #print 'herbivore feeding... gain = ', lost_biomass_prey
            if omnivore:
                biomass_gained = lost_biomass_prey*OMNIVORY_TRADEOFF #*HERBIVORY_EFFICIENCY
            else:
                if self.mutualist:
                    biomass_gained = lost_biomass_prey*MUTUALISTIC_FRACTION
                else:
                    biomass_gained = lost_biomass_prey*HERBIVORY_EFFICIENCY
            #if not self.mutualist:
            prey.resource -= lost_biomass_prey
        else:
#            print 'predation'
#            if prey.mutualist:
#                print 'a mutualist has been eaten'
            if self.mutualist:
                print 'I am a mutualist predator', self.species_id
            
            #print 'predation'
            
            biomass_gained = prey.resource*EFFICIENCY_TRANSFER
            
        self.resource += biomass_gained
        
        if self.resource > MAX_RESOURCE:
            self.resource = MAX_RESOURCE
        
      
    def ready_to_mate(self):
        if self.resource > MAX_RESOURCE*MATING_RESOURCE:
            return True
        
        return False
    
    def reproduce(self):
        resource_given = self.resource*MATING_ENERGY      
        self.resource -= resource_given
        
#        if self.mutualist:
#            print 'mutualist reproducing', self.species_id, 'resource lost', resource_given
#        
        return resource_given
    
    def synthesis(self):
        synthesised_resource = self.resource*SYNTHESIS_ABILITY
        self.resource += synthesised_resource
        
        if self.mutualist:
            print 'I am an autotroph mutualist'
            print 'species id', self.species_id
        
        if self.resource > MAX_RESOURCE:
            self.resource = MAX_RESOURCE
        
    def set_mutualistic_partner(self, partner):
        if self.mutualist == False:
            return
        self.current_host = partner.species_id
        self.mut_efficiency = MUTUALISTIC_EFFICIENCY
    
    def mutualistic_cool_off(self):
        self.mut_efficiency = self.mut_efficiency*MUTUALISTIC_COOLING 
        if self.mut_efficiency < MIN_MUTUALISTIC_EFF:
            self.current_host = None
            
    def reset_mutualistic_state(self):
        self.current_host = None
    
    
        
    