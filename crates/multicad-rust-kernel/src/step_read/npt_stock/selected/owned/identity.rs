//! Exact defining STEP mapping leaves and incidence, captured before AST drop.
use super::*;
struct Sink<'a> {
    length: usize,
    bytes: Option<&'a mut Vec<u8>>,
}
impl Sink<'_> {
    fn put(&mut self, bytes: &[u8], budget: &mut Budget<'_>) -> Result<()> {
        budget.step(bytes.len() + 8)?;
        self.length = self.length.checked_add(bytes.len()).ok_or(Error::Limit)?;
        if let Some(output) = &mut self.bytes {
            output.extend_from_slice(bytes);
        }
        Ok(())
    }
    fn number(&mut self, value: usize, budget: &mut Budget<'_>) -> Result<()> {
        self.put(&(value as u64).to_le_bytes(), budget)
    }
    fn text(&mut self, value: &str, budget: &mut Budget<'_>) -> Result<()> {
        self.number(value.len(), budget)?;
        self.put(value.as_bytes(), budget)
    }
}
fn reference(id: u32, nodes: &mut Vec<u32>, budget: &mut Budget<'_>) -> Result<usize> {
    budget.step(32)?;
    if let Some(index) = nodes.iter().position(|v| *v == id) {
        return Ok(index);
    }
    // The admitted Cartesian operator and source AXIS2 frame contain only
    // their finite point/direction records. Keep metadata traversal bounded
    // independently of the unrelated represented body and its large tensors.
    if nodes.len() == 32 {
        return Err(Error::Limit);
    }
    let index = nodes.len();
    nodes.push(id);
    Ok(index)
}
fn value(
    input: &Value<'_>,
    numbers: &mut numeric::Numbers<'_>,
    nodes: &mut Vec<u32>,
    sink: &mut Sink<'_>,
    budget: &mut Budget<'_>,
) -> Result<()> {
    match input {
        // Labels are not defining mapping geometry. Their omission also
        // prevents an innocuous display-name edit from rebinding topology.
        Value::String(_) => sink.put(b"s", budget),
        Value::Number(_) => {
            sink.put(b"n", budget)?;
            sink.text(numbers.scalar(input, budget)?.spelling(), budget)
        }
        Value::Reference(id) => {
            sink.put(b"r", budget)?;
            sink.number(reference(*id, nodes, budget)?, budget)
        }
        Value::List(items) => {
            sink.put(b"l", budget)?;
            sink.number(items.len(), budget)?;
            for item in items {
                value(item, numbers, nodes, sink, budget)?;
            }
            Ok(())
        }
        Value::Omitted => sink.put(b"$", budget),
        Value::Derived => sink.put(b"*", budget),
        Value::Enum(name) => {
            sink.put(b"e", budget)?;
            sink.text(name, budget)
        }
        Value::Typed(_) => Err(Error::Unsupported),
    }
}
fn encode(
    arena: &BTreeMap<u32, Vec<Entity<'_>>>,
    nodes: &mut Vec<u32>,
    sink: &mut Sink<'_>,
    budget: &mut Budget<'_>,
) -> Result<()> {
    let mut index = 0;
    while index < nodes.len() {
        budget.step(32)?;
        let components = arena.get(&nodes[index]).ok_or(Error::Invalid)?;
        sink.number(index, budget)?;
        sink.number(components.len(), budget)?;
        for entity in components {
            sink.text(entity.name, budget)?;
            sink.number(entity.args.len(), budget)?;
            let mut numbers = numeric::Numbers::new(entity, budget)?;
            for arg in &entity.args {
                value(arg, &mut numbers, nodes, sink, budget)?;
            }
            numbers.finish()?;
        }
        index += 1;
    }
    sink.number(nodes.len(), budget)
}
pub(super) fn capture(selected: &Selected<'_, '_>, budget: &mut Budget<'_>) -> Result<Vec<u8>> {
    let Some(link) = selected.member.mapping_link() else {
        return budget.allocate(0);
    };
    let mut nodes = budget.allocate(32)?;
    reference(link.operator_entity(), &mut nodes, budget)?;
    reference(link.source_origin_entity(), &mut nodes, budget)?;
    let mut count = Sink {
        length: 0,
        bytes: None,
    };
    encode(
        selected.parsed.document.records(),
        &mut nodes,
        &mut count,
        budget,
    )?;
    let mut output = budget.allocate(count.length)?;
    let mut sink = Sink {
        length: 0,
        bytes: Some(&mut output),
    };
    encode(
        selected.parsed.document.records(),
        &mut nodes,
        &mut sink,
        budget,
    )?;
    assert_eq!(
        count.length, sink.length,
        "immutable mapping graph encoding"
    );
    budget.release(&nodes);
    Ok(output)
}
