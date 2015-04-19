$(document).on('click', '.selCell tbody td', function(){
	var el = $(this);
	$(this).toggleClass("cellsSelected", this.clicked);
	el.trigger("change");
});	
var selectCellBinding = new Shiny.InputBinding();
$.extend(selectCellBinding, {
	find: function(scope) {
		return $(scope).find(".selCell");;
	},
	getValue: function(el){
		var rowIndex = $(el).children().children().children('.cellsSelected').map(function() { var ri = $(this).parent().index();
		if (ri ==-1)
			return -1;
		return ri;}).get();
		var colIndex = $(el).children().children().children('.cellsSelected').map(function() {
                 return $(this).index();
              }).get();
		return [[rowIndex], [colIndex]];
	},
	setValue: function(el, value) {
	},
	subscribe: function(el, callback) {
		$(el).on("change.selectCellBinding", function(e) {
			callback();
		});
	},
	unsubscribe: function(el) {
	  $(el).off(".selectCellBinding");
	}
});
Shiny.inputBindings.register(selectCellBinding);
